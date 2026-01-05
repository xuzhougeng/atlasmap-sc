// Package service provides business logic for the tile server.
package service

import (
	"context"
	"fmt"
	"math"
	"math/rand"
	"sort"

	"github.com/soma-tiles/server/internal/data/soma"
	"github.com/soma-tiles/server/internal/destore"
)

// DEService handles differential expression analysis.
type DEService struct {
	registry interface {
		Get(datasetID string) *TileService
	}
}

// NewDEService creates a new DE service.
func NewDEService(registry interface{ Get(datasetID string) *TileService }) *DEService {
	return &DEService{registry: registry}
}

// ExecuteDEJob runs the DE analysis for a job (called by JobManager worker).
func (s *DEService) ExecuteDEJob(ctx context.Context, store *destore.Store, jobID string) error {
	// Load job from store
	job, err := store.GetJob(jobID)
	if err != nil {
		return fmt.Errorf("failed to get job: %w", err)
	}
	if job == nil {
		return fmt.Errorf("job not found: %s", jobID)
	}

	svc := s.registry.Get(job.Params.DatasetID)
	if svc == nil {
		return fmt.Errorf("dataset not found: %s", job.Params.DatasetID)
	}
	sr := svc.Soma()
	if sr == nil {
		return fmt.Errorf("soma not configured for dataset: %s", job.Params.DatasetID)
	}
	if !sr.Supported() {
		return soma.ErrUnsupported
	}

	// Phase 1: Get cell groups
	store.UpdateJobProgress(jobID, "loading_obs", 0, 100)

	groupIdx, err := sr.ObsGroupIndex(job.Params.Groupby)
	if err != nil {
		return fmt.Errorf("failed to load obs groupby: %w", err)
	}

	// Collect cells for group1
	var group1Cells []int64
	for _, g := range job.Params.Group1 {
		if cells, ok := groupIdx[g]; ok {
			group1Cells = append(group1Cells, cells...)
		}
	}
	if len(group1Cells) == 0 {
		return fmt.Errorf("group1 has no cells for values: %v", job.Params.Group1)
	}

	// Collect cells for group2 (if not specified, use "rest")
	var group2Cells []int64
	if len(job.Params.Group2) == 0 {
		// One-vs-rest: all cells not in group1
		g1Set := make(map[int64]bool, len(group1Cells))
		for _, c := range group1Cells {
			g1Set[c] = true
		}
		for _, cells := range groupIdx {
			for _, c := range cells {
				if !g1Set[c] {
					group2Cells = append(group2Cells, c)
				}
			}
		}
	} else {
		for _, g := range job.Params.Group2 {
			if cells, ok := groupIdx[g]; ok {
				group2Cells = append(group2Cells, cells...)
			}
		}
	}
	if len(group2Cells) == 0 {
		return fmt.Errorf("group2 has no cells")
	}

	// Sample cells
	rng := rand.New(rand.NewSource(int64(job.Params.Seed)))
	n1 := len(group1Cells)
	n2 := len(group2Cells)
	maxCells := job.Params.MaxCellsPerGroup
	if n1 > maxCells {
		rng.Shuffle(len(group1Cells), func(i, j int) {
			group1Cells[i], group1Cells[j] = group1Cells[j], group1Cells[i]
		})
		group1Cells = group1Cells[:maxCells]
		n1 = maxCells
	}
	if n2 > maxCells {
		rng.Shuffle(len(group2Cells), func(i, j int) {
			group2Cells[i], group2Cells[j] = group2Cells[j], group2Cells[i]
		})
		group2Cells = group2Cells[:maxCells]
		n2 = maxCells
	}

	// Update n1/n2 in DB
	store.UpdateJobCounts(jobID, n1, n2)

	if ctx.Err() != nil {
		return ctx.Err()
	}

	// Phase 2: Get gene map
	store.UpdateJobProgress(jobID, "loading_genes", 0, 100)
	geneMap, err := sr.AllGenes()
	if err != nil {
		return fmt.Errorf("failed to load gene map: %w", err)
	}
	// Invert map: joinid -> gene_id
	geneID := make(map[int64]string, len(geneMap))
	for g, id := range geneMap {
		geneID[id] = g
	}
	nGenes := len(geneMap)

	// Phase 3: Scan X for group1
	store.UpdateJobProgress(jobID, "scanning_group1", 0, nGenes)

	type geneAcc struct {
		sum   float64
		sumsq float64
		nnz   int
		vals  []float32
	}
	acc1 := make(map[int64]*geneAcc)
	err = sr.ScanXForCells(group1Cells, func(cell, gene int64, val float32) {
		a, ok := acc1[gene]
		if !ok {
			a = &geneAcc{}
			acc1[gene] = a
		}
		v := float64(val)
		a.sum += v
		a.sumsq += v * v
		a.nnz++
		a.vals = append(a.vals, val)
	})
	if err != nil {
		return fmt.Errorf("failed to scan X for group1: %w", err)
	}

	if ctx.Err() != nil {
		return ctx.Err()
	}

	// Phase 4: Scan X for group2
	store.UpdateJobProgress(jobID, "scanning_group2", 0, nGenes)

	acc2 := make(map[int64]*geneAcc)
	err = sr.ScanXForCells(group2Cells, func(cell, gene int64, val float32) {
		a, ok := acc2[gene]
		if !ok {
			a = &geneAcc{}
			acc2[gene] = a
		}
		v := float64(val)
		a.sum += v
		a.sumsq += v * v
		a.nnz++
		a.vals = append(a.vals, val)
	})
	if err != nil {
		return fmt.Errorf("failed to scan X for group2: %w", err)
	}

	if ctx.Err() != nil {
		return ctx.Err()
	}

	// Phase 5: Compute statistics
	store.UpdateJobProgress(jobID, "computing_stats", 0, nGenes)

	wantTtest := contains(job.Params.Tests, "ttest")
	wantRanksum := contains(job.Params.Tests, "ranksum")

	type geneResult struct {
		geneJoinID int64
		gene       string
		mean1      float64
		mean2      float64
		pct1       float64
		pct2       float64
		log2fc     float64
		pTtest     float64
		pRanksum   float64
	}

	results := make([]geneResult, 0, nGenes)

	for gid, gname := range geneID {
		a1 := acc1[gid]
		a2 := acc2[gid]

		mean1, mean2 := 0.0, 0.0
		var1, var2 := 0.0, 0.0
		pct1, pct2 := 0.0, 0.0
		nnz1, nnz2 := 0, 0

		if a1 != nil && a1.nnz > 0 {
			mean1 = a1.sum / float64(n1)
			pct1 = float64(a1.nnz) / float64(n1)
			nnz1 = a1.nnz
			if n1 > 1 {
				sumSq := a1.sumsq
				var1 = (sumSq/float64(n1) - mean1*mean1) * float64(n1) / float64(n1-1)
			}
		}
		if a2 != nil && a2.nnz > 0 {
			mean2 = a2.sum / float64(n2)
			pct2 = float64(a2.nnz) / float64(n2)
			nnz2 = a2.nnz
			if n2 > 1 {
				sumSq := a2.sumsq
				var2 = (sumSq/float64(n2) - mean2*mean2) * float64(n2) / float64(n2-1)
			}
		}

		log2fc := 0.0
		eps := 1e-9
		if mean2 > eps || mean1 > eps {
			log2fc = math.Log2((mean1 + eps) / (mean2 + eps))
		}

		pTtest := 1.0
		if wantTtest && (nnz1 > 0 || nnz2 > 0) {
			pTtest = welchTTest(mean1, var1, n1, mean2, var2, n2)
		}

		pRanksum := 1.0
		if wantRanksum && (nnz1 > 0 || nnz2 > 0) {
			var vals1, vals2 []float32
			if a1 != nil {
				vals1 = a1.vals
			}
			if a2 != nil {
				vals2 = a2.vals
			}
			pRanksum = mannWhitneyU(vals1, n1, vals2, n2)
		}

		results = append(results, geneResult{
			geneJoinID: gid,
			gene:       gname,
			mean1:      mean1,
			mean2:      mean2,
			pct1:       pct1,
			pct2:       pct2,
			log2fc:     log2fc,
			pTtest:     pTtest,
			pRanksum:   pRanksum,
		})
	}

	if ctx.Err() != nil {
		return ctx.Err()
	}

	// Phase 6: FDR correction
	store.UpdateJobProgress(jobID, "computing_fdr", 0, nGenes)

	pTtests := make([]float64, len(results))
	pRanksums := make([]float64, len(results))
	for i, r := range results {
		pTtests[i] = r.pTtest
		pRanksums[i] = r.pRanksum
	}

	fdrTtests := benjaminiHochberg(pTtests)
	fdrRanksums := benjaminiHochberg(pRanksums)

	// Build final results for DB
	items := make([]*destore.DEGeneResult, len(results))
	for i, r := range results {
		items[i] = &destore.DEGeneResult{
			Gene:       r.gene,
			GeneJoinID: r.geneJoinID,
			Mean1:      r.mean1,
			Mean2:      r.mean2,
			Pct1:       r.pct1,
			Pct2:       r.pct2,
			Log2FC:     r.log2fc,
			PTtest:     r.pTtest,
			FDRTtest:   fdrTtests[i],
			PRanksum:   r.pRanksum,
			FDRRanksum: fdrRanksums[i],
		}
	}

	// Sort by FDR before insert (for consistent ordering)
	sort.Slice(items, func(i, j int) bool {
		if items[i].FDRRanksum != items[j].FDRRanksum {
			return items[i].FDRRanksum < items[j].FDRRanksum
		}
		return math.Abs(items[i].Log2FC) > math.Abs(items[j].Log2FC)
	})

	// Phase 7: Write results to DB
	store.UpdateJobProgress(jobID, "saving_results", 0, len(items))

	if err := store.InsertResults(jobID, items); err != nil {
		return fmt.Errorf("failed to save results: %w", err)
	}

	return nil
}

// welchTTest computes the p-value for Welch's t-test (two-tailed).
func welchTTest(mean1, var1 float64, n1 int, mean2, var2 float64, n2 int) float64 {
	if n1 < 2 || n2 < 2 {
		return 1.0
	}
	if var1 <= 0 && var2 <= 0 {
		if mean1 == mean2 {
			return 1.0
		}
		return 0.0
	}

	se1 := var1 / float64(n1)
	se2 := var2 / float64(n2)
	seDiff := math.Sqrt(se1 + se2)
	if seDiff < 1e-15 {
		if mean1 == mean2 {
			return 1.0
		}
		return 0.0
	}

	t := (mean1 - mean2) / seDiff

	num := (se1 + se2) * (se1 + se2)
	den := 0.0
	if n1 > 1 && se1 > 0 {
		den += se1 * se1 / float64(n1-1)
	}
	if n2 > 1 && se2 > 0 {
		den += se2 * se2 / float64(n2-1)
	}
	if den < 1e-15 {
		return 1.0
	}
	df := num / den
	if df < 1 {
		df = 1
	}

	return 2 * studentTCDF(-math.Abs(t), df)
}

func studentTCDF(t, df float64) float64 {
	if df <= 0 {
		return 0.5
	}
	x := df / (df + t*t)
	beta := incompleteBeta(x, df/2, 0.5)
	if t < 0 {
		return 0.5 * beta
	}
	return 1 - 0.5*beta
}

func incompleteBeta(x, a, b float64) float64 {
	if x <= 0 {
		return 0
	}
	if x >= 1 {
		return 1
	}

	if x > (a+1)/(a+b+2) {
		return 1 - incompleteBeta(1-x, b, a)
	}

	const maxIter = 200
	const eps = 1e-10

	lnBeta := lgamma(a) + lgamma(b) - lgamma(a+b)
	front := math.Exp(a*math.Log(x) + b*math.Log(1-x) - lnBeta)

	f := 1.0
	c := 1.0
	d := 0.0

	for i := 0; i <= maxIter; i++ {
		m := float64(i / 2)
		var num float64
		if i == 0 {
			num = 1.0
		} else if i%2 == 0 {
			num = m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m))
		} else {
			num = -((a + m) * (a + b + m) * x) / ((a + 2*m) * (a + 2*m + 1))
		}

		d = 1 + num*d
		if math.Abs(d) < eps {
			d = eps
		}
		d = 1 / d

		c = 1 + num/c
		if math.Abs(c) < eps {
			c = eps
		}

		f *= d * c
		if math.Abs(d*c-1) < eps {
			break
		}
	}

	return front * (f - 1) / a
}

func lgamma(x float64) float64 {
	g, _ := math.Lgamma(x)
	return g
}

func mannWhitneyU(vals1 []float32, n1 int, vals2 []float32, n2 int) float64 {
	if n1 == 0 || n2 == 0 {
		return 1.0
	}

	type entry struct {
		val   float32
		group int
	}
	combined := make([]entry, 0, len(vals1)+len(vals2)+(n1-len(vals1))+(n2-len(vals2)))

	for _, v := range vals1 {
		combined = append(combined, entry{val: v, group: 1})
	}
	for _, v := range vals2 {
		combined = append(combined, entry{val: v, group: 2})
	}
	nZeros1 := n1 - len(vals1)
	nZeros2 := n2 - len(vals2)
	for i := 0; i < nZeros1; i++ {
		combined = append(combined, entry{val: 0, group: 1})
	}
	for i := 0; i < nZeros2; i++ {
		combined = append(combined, entry{val: 0, group: 2})
	}

	sort.Slice(combined, func(i, j int) bool {
		return combined[i].val < combined[j].val
	})

	N := len(combined)
	ranks := make([]float64, N)
	i := 0
	for i < N {
		j := i
		for j < N && combined[j].val == combined[i].val {
			j++
		}
		avgRank := float64(i+j+1) / 2.0
		for k := i; k < j; k++ {
			ranks[k] = avgRank
		}
		i = j
	}

	R1 := 0.0
	for i, e := range combined {
		if e.group == 1 {
			R1 += ranks[i]
		}
	}

	n1f := float64(n1)
	n2f := float64(n2)
	U1 := R1 - n1f*(n1f+1)/2
	U2 := n1f*n2f - U1
	U := math.Min(U1, U2)

	muU := n1f * n2f / 2

	tieSum := 0.0
	i = 0
	for i < N {
		j := i
		for j < N && combined[j].val == combined[i].val {
			j++
		}
		t := float64(j - i)
		if t > 1 {
			tieSum += t*t*t - t
		}
		i = j
	}

	Nf := float64(N)
	sigmaU := math.Sqrt(n1f * n2f * ((Nf + 1) - tieSum/(Nf*(Nf-1))) / 12)

	if sigmaU < 1e-10 {
		return 1.0
	}

	z := (U - muU + 0.5) / sigmaU
	p := 2 * normalCDF(-math.Abs(z))
	return p
}

func normalCDF(x float64) float64 {
	return 0.5 * (1 + math.Erf(x/math.Sqrt2))
}

func benjaminiHochberg(pvals []float64) []float64 {
	n := len(pvals)
	if n == 0 {
		return nil
	}

	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(i, j int) bool {
		return pvals[idx[i]] < pvals[idx[j]]
	})

	fdr := make([]float64, n)
	minP := 1.0
	for i := n - 1; i >= 0; i-- {
		origIdx := idx[i]
		rank := i + 1
		adjusted := pvals[origIdx] * float64(n) / float64(rank)
		if adjusted > 1 {
			adjusted = 1
		}
		if adjusted < minP {
			minP = adjusted
		} else {
			adjusted = minP
		}
		fdr[origIdx] = adjusted
	}

	return fdr
}

func contains(slice []string, item string) bool {
	for _, s := range slice {
		if s == item {
			return true
		}
	}
	return false
}
