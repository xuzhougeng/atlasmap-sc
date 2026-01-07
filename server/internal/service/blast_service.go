// Package service provides business logic for the tile server.
package service

import (
	"bufio"
	"context"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/soma-tiles/server/internal/blaststore"
)

// BlastServiceRegistry provides access to dataset metadata for gene validation.
type BlastServiceRegistry interface {
	DatasetIDs() []string
	BlastPPath(datasetID string) string
	Get(datasetID string) *TileService
}

// BlastService handles BLASTP searches across datasets.
type BlastService struct {
	registry BlastServiceRegistry
}

// NewBlastService creates a new BLAST service.
func NewBlastService(registry BlastServiceRegistry) *BlastService {
	return &BlastService{registry: registry}
}

// ExecuteBlastJob runs the BLAST search for a job (called by BlastJobManager worker).
func (s *BlastService) ExecuteBlastJob(ctx context.Context, store *blaststore.Store, jobID string) error {
	// Load job from store
	job, err := store.GetJob(jobID)
	if err != nil {
		return fmt.Errorf("failed to get job: %w", err)
	}
	if job == nil {
		return fmt.Errorf("job not found: %s", jobID)
	}

	// Collect all datasets with blastp_path
	type dbInfo struct {
		path     string
		datasets []string
	}
	dbMap := make(map[string]*dbInfo) // blastp_path -> info

	targetDatasets := job.Params.Datasets
	if len(targetDatasets) == 0 {
		// Use all datasets
		targetDatasets = s.registry.DatasetIDs()
	}

	for _, dsID := range targetDatasets {
		path := s.registry.BlastPPath(dsID)
		if path == "" {
			continue
		}
		if info, ok := dbMap[path]; ok {
			info.datasets = append(info.datasets, dsID)
		} else {
			dbMap[path] = &dbInfo{path: path, datasets: []string{dsID}}
		}
	}

	if len(dbMap) == 0 {
		return fmt.Errorf("no datasets have blastp_path configured")
	}

	// Apply defaults
	maxHits := job.Params.MaxHits
	if maxHits <= 0 {
		maxHits = 10
	}
	evalue := job.Params.Evalue
	if evalue <= 0 {
		evalue = 1e-5
	}
	numThreads := job.Params.NumThreads
	if numThreads <= 0 {
		numThreads = 1
	}

	// Create temp file for query sequence
	tmpFile, err := os.CreateTemp("", "blast_query_*.fasta")
	if err != nil {
		return fmt.Errorf("failed to create temp file: %w", err)
	}
	tmpPath := tmpFile.Name()
	defer os.Remove(tmpPath)

	// Write query sequence (wrap in FASTA format if needed)
	seq := strings.TrimSpace(job.Params.Sequence)
	if !strings.HasPrefix(seq, ">") {
		seq = ">query\n" + seq
	}
	if _, err := tmpFile.WriteString(seq); err != nil {
		tmpFile.Close()
		return fmt.Errorf("failed to write query sequence: %w", err)
	}
	tmpFile.Close()

	// Update progress
	store.UpdateJobProgress(jobID, "running_blast", 0, len(dbMap))

	// Run BLAST for each unique database in parallel (with limit)
	type blastResult struct {
		dbPath string
		hits   []rawBlastHit
		err    error
	}

	resultCh := make(chan blastResult, len(dbMap))
	var wg sync.WaitGroup

	// Limit concurrency to avoid spawning too many blastp processes
	sem := make(chan struct{}, 4) // max 4 concurrent blastp

	for dbPath := range dbMap {
		wg.Add(1)
		go func(dbPath string) {
			defer wg.Done()
			sem <- struct{}{}
			defer func() { <-sem }()

			hits, err := runBlastP(ctx, dbPath, tmpPath, maxHits, evalue, numThreads)
			resultCh <- blastResult{dbPath: dbPath, hits: hits, err: err}
		}(dbPath)
	}

	// Wait for all to complete
	go func() {
		wg.Wait()
		close(resultCh)
	}()

	// Collect results
	var allHits []*blaststore.BlastHit
	done := 0
	for res := range resultCh {
		done++
		store.UpdateJobProgress(jobID, "running_blast", done, len(dbMap))

		if ctx.Err() != nil {
			return ctx.Err()
		}

		if res.err != nil {
			// Log error but continue with other databases
			continue
		}

		info := dbMap[res.dbPath]
		for _, hit := range res.hits {
			// Expand hit to all datasets using this database
			for _, dsID := range info.datasets {
				// Optional: validate gene exists in dataset
				svc := s.registry.Get(dsID)
				if svc != nil {
					md := svc.Metadata()
					if md.GeneIndex != nil {
						if _, ok := md.GeneIndex[hit.sseqid]; !ok {
							// Gene not in this dataset's gene index, skip
							continue
						}
					}
				}

				allHits = append(allHits, &blaststore.BlastHit{
					DatasetID: dsID,
					GeneID:    hit.sseqid,
					Pident:    hit.pident,
					Length:    hit.length,
					Evalue:    hit.evalue,
					Bitscore:  hit.bitscore,
				})
			}
		}
	}

	if ctx.Err() != nil {
		return ctx.Err()
	}

	// Store results
	store.UpdateJobProgress(jobID, "saving_results", 0, len(allHits))
	if err := store.InsertResults(jobID, allHits); err != nil {
		return fmt.Errorf("failed to save results: %w", err)
	}

	return nil
}

// rawBlastHit represents a single hit from blastp output
type rawBlastHit struct {
	sseqid   string
	pident   float64
	length   int
	evalue   float64
	bitscore float64
}

// runBlastP runs blastp and parses the output
func runBlastP(ctx context.Context, dbPath, queryPath string, maxHits int, evalue float64, numThreads int) ([]rawBlastHit, error) {
	// Check if blastp exists
	blastpPath, err := exec.LookPath("blastp")
	if err != nil {
		return nil, fmt.Errorf("blastp not found in PATH: %w", err)
	}

	// Verify database files exist
	dbDir := filepath.Dir(dbPath)
	dbBase := filepath.Base(dbPath)
	found := false
	for _, ext := range []string{".phr", ".pin", ".psq", ".pal"} {
		if _, err := os.Stat(dbPath + ext); err == nil {
			found = true
			break
		}
	}
	if !found {
		return nil, fmt.Errorf("BLAST database not found at %s (checked in %s for %s.*)", dbPath, dbDir, dbBase)
	}

	// Build blastp command
	// outfmt 6: tab-separated with fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
	// We only need: sseqid pident length evalue bitscore
	args := []string{
		"-db", dbPath,
		"-query", queryPath,
		"-outfmt", "6 sseqid pident length evalue bitscore",
		"-max_target_seqs", strconv.Itoa(maxHits),
		"-evalue", strconv.FormatFloat(evalue, 'e', 2, 64),
		"-num_threads", strconv.Itoa(numThreads),
	}

	cmd := exec.CommandContext(ctx, blastpPath, args...)
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return nil, fmt.Errorf("failed to create stdout pipe: %w", err)
	}

	if err := cmd.Start(); err != nil {
		return nil, fmt.Errorf("failed to start blastp: %w", err)
	}

	// Parse output
	var hits []rawBlastHit
	scanner := bufio.NewScanner(stdout)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		pident, _ := strconv.ParseFloat(fields[1], 64)
		length, _ := strconv.Atoi(fields[2])
		ev, _ := strconv.ParseFloat(fields[3], 64)
		bitscore, _ := strconv.ParseFloat(fields[4], 64)

		hits = append(hits, rawBlastHit{
			sseqid:   fields[0],
			pident:   pident,
			length:   length,
			evalue:   ev,
			bitscore: bitscore,
		})
	}

	if err := cmd.Wait(); err != nil {
		// Check if cancelled
		if ctx.Err() != nil {
			return nil, ctx.Err()
		}
		// blastp returns non-zero for no hits, which is fine
		// Only fail on actual errors
		if exitErr, ok := err.(*exec.ExitError); ok && exitErr.ExitCode() != 0 {
			// Return empty hits rather than error for "no hits"
			return hits, nil
		}
		return nil, fmt.Errorf("blastp failed: %w", err)
	}

	return hits, nil
}

