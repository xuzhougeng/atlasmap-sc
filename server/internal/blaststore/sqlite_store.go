// Package blaststore provides persistent storage for BLASTP job state and results using SQLite.
package blaststore

import (
	"database/sql"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"time"

	_ "modernc.org/sqlite"
)

// JobStatus represents the current state of a BLAST job.
type JobStatus string

const (
	JobStatusQueued    JobStatus = "queued"
	JobStatusRunning   JobStatus = "running"
	JobStatusCompleted JobStatus = "completed"
	JobStatusFailed    JobStatus = "failed"
	JobStatusCancelled JobStatus = "cancelled"
)

// BlastJobParams contains the parameters for a BLAST job.
type BlastJobParams struct {
	Sequence    string   `json:"sequence"`     // Query protein sequence (FASTA or raw)
	MaxHits     int      `json:"max_hits"`     // Max hits per database (default 10)
	Evalue      float64  `json:"evalue"`       // E-value cutoff (default 1e-5)
	Datasets    []string `json:"datasets"`     // Optional: limit to specific datasets (empty = all)
	NumThreads  int      `json:"num_threads"`  // Threads per blastp process (default 1)
}

// BlastJobProgress represents the progress of a BLAST job.
type BlastJobProgress struct {
	Phase string `json:"phase"`
	Done  int    `json:"done"`
	Total int    `json:"total"`
}

// BlastJob represents a BLASTP search job.
type BlastJob struct {
	ID         string           `json:"job_id"`
	Status     JobStatus        `json:"status"`
	Params     BlastJobParams   `json:"params"`
	Progress   BlastJobProgress `json:"progress"`
	CreatedAt  time.Time        `json:"created_at"`
	StartedAt  *time.Time       `json:"started_at,omitempty"`
	FinishedAt *time.Time       `json:"finished_at,omitempty"`
	Error      string           `json:"error,omitempty"`
}

// BlastHit contains a single BLASTP hit result.
type BlastHit struct {
	DatasetID string  `json:"dataset_id"`
	GeneID    string  `json:"gene_id"`
	Pident    float64 `json:"pident"`   // Percent identity
	Length    int     `json:"length"`   // Alignment length
	Evalue    float64 `json:"evalue"`   // E-value
	Bitscore  float64 `json:"bitscore"` // Bit score
}

// Store provides persistent storage for BLAST jobs using SQLite.
type Store struct {
	db *sql.DB
	mu sync.Mutex
}

// NewStore creates a new SQLite-based BLAST store.
func NewStore(dbPath string) (*Store, error) {
	// Ensure directory exists
	dir := filepath.Dir(dbPath)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return nil, fmt.Errorf("failed to create directory for sqlite: %w", err)
	}

	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open sqlite: %w", err)
	}

	// Enable WAL mode for better concurrency
	if _, err := db.Exec("PRAGMA journal_mode=WAL"); err != nil {
		db.Close()
		return nil, fmt.Errorf("failed to enable WAL: %w", err)
	}

	s := &Store{db: db}
	if err := s.migrate(); err != nil {
		db.Close()
		return nil, fmt.Errorf("failed to migrate: %w", err)
	}

	return s, nil
}

// Close closes the database connection.
func (s *Store) Close() error {
	return s.db.Close()
}

func (s *Store) migrate() error {
	schema := `
	CREATE TABLE IF NOT EXISTS blast_jobs (
		job_id TEXT PRIMARY KEY,
		status TEXT NOT NULL,
		params_json TEXT NOT NULL,
		phase TEXT DEFAULT '',
		done INTEGER DEFAULT 0,
		total INTEGER DEFAULT 0,
		error TEXT DEFAULT '',
		created_at TEXT NOT NULL,
		started_at TEXT,
		finished_at TEXT
	);

	CREATE INDEX IF NOT EXISTS idx_blast_jobs_status ON blast_jobs(status);
	CREATE INDEX IF NOT EXISTS idx_blast_jobs_finished ON blast_jobs(finished_at);

	CREATE TABLE IF NOT EXISTS blast_results (
		id INTEGER PRIMARY KEY AUTOINCREMENT,
		job_id TEXT NOT NULL,
		dataset_id TEXT NOT NULL,
		gene_id TEXT NOT NULL,
		pident REAL NOT NULL,
		length INTEGER NOT NULL,
		evalue REAL NOT NULL,
		bitscore REAL NOT NULL,
		FOREIGN KEY (job_id) REFERENCES blast_jobs(job_id) ON DELETE CASCADE
	);

	CREATE INDEX IF NOT EXISTS idx_blast_results_job ON blast_results(job_id);
	CREATE INDEX IF NOT EXISTS idx_blast_results_job_evalue ON blast_results(job_id, evalue);
	CREATE INDEX IF NOT EXISTS idx_blast_results_job_bitscore ON blast_results(job_id, bitscore DESC);
	`
	_, err := s.db.Exec(schema)
	return err
}

// CreateJob creates a new job record with status=queued.
func (s *Store) CreateJob(job *BlastJob) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	paramsJSON, err := json.Marshal(job.Params)
	if err != nil {
		return fmt.Errorf("failed to marshal params: %w", err)
	}

	_, err = s.db.Exec(`
		INSERT INTO blast_jobs (job_id, status, params_json, phase, done, total, error, created_at, started_at, finished_at)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	`,
		job.ID,
		string(job.Status),
		string(paramsJSON),
		job.Progress.Phase,
		job.Progress.Done,
		job.Progress.Total,
		job.Error,
		job.CreatedAt.Format(time.RFC3339),
		nil,
		nil,
	)
	return err
}

// GetJob retrieves a job by ID.
func (s *Store) GetJob(jobID string) (*BlastJob, error) {
	row := s.db.QueryRow(`
		SELECT job_id, status, params_json, phase, done, total, error, created_at, started_at, finished_at
		FROM blast_jobs WHERE job_id = ?
	`, jobID)

	var job BlastJob
	var paramsJSON string
	var createdAtStr string
	var startedAtStr, finishedAtStr sql.NullString

	err := row.Scan(
		&job.ID,
		&job.Status,
		&paramsJSON,
		&job.Progress.Phase,
		&job.Progress.Done,
		&job.Progress.Total,
		&job.Error,
		&createdAtStr,
		&startedAtStr,
		&finishedAtStr,
	)
	if err == sql.ErrNoRows {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}

	if err := json.Unmarshal([]byte(paramsJSON), &job.Params); err != nil {
		return nil, fmt.Errorf("failed to unmarshal params: %w", err)
	}

	job.CreatedAt, _ = time.Parse(time.RFC3339, createdAtStr)
	if startedAtStr.Valid {
		t, _ := time.Parse(time.RFC3339, startedAtStr.String)
		job.StartedAt = &t
	}
	if finishedAtStr.Valid {
		t, _ := time.Parse(time.RFC3339, finishedAtStr.String)
		job.FinishedAt = &t
	}

	return &job, nil
}

// UpdateJobStatus updates the job status and optional fields.
func (s *Store) UpdateJobStatus(jobID string, status JobStatus, errMsg string) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	var finishedAt *string
	if status == JobStatusCompleted || status == JobStatusFailed || status == JobStatusCancelled {
		t := time.Now().Format(time.RFC3339)
		finishedAt = &t
	}

	_, err := s.db.Exec(`
		UPDATE blast_jobs SET status = ?, error = ?, finished_at = COALESCE(?, finished_at)
		WHERE job_id = ?
	`, string(status), errMsg, finishedAt, jobID)
	return err
}

// UpdateJobStarted marks a job as running with start time.
func (s *Store) UpdateJobStarted(jobID string) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	now := time.Now().Format(time.RFC3339)
	_, err := s.db.Exec(`
		UPDATE blast_jobs SET status = ?, started_at = ?
		WHERE job_id = ?
	`, string(JobStatusRunning), now, jobID)
	return err
}

// UpdateJobProgress updates the progress fields.
func (s *Store) UpdateJobProgress(jobID string, phase string, done, total int) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	_, err := s.db.Exec(`
		UPDATE blast_jobs SET phase = ?, done = ?, total = ?
		WHERE job_id = ?
	`, phase, done, total, jobID)
	return err
}

// InsertResults inserts BLAST hit results in a batch transaction.
func (s *Store) InsertResults(jobID string, results []*BlastHit) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	tx, err := s.db.Begin()
	if err != nil {
		return err
	}
	defer tx.Rollback()

	stmt, err := tx.Prepare(`
		INSERT INTO blast_results (job_id, dataset_id, gene_id, pident, length, evalue, bitscore)
		VALUES (?, ?, ?, ?, ?, ?, ?)
	`)
	if err != nil {
		return err
	}
	defer stmt.Close()

	for _, r := range results {
		_, err := stmt.Exec(
			jobID, r.DatasetID, r.GeneID,
			r.Pident, r.Length, r.Evalue, r.Bitscore,
		)
		if err != nil {
			return err
		}
	}

	return tx.Commit()
}

// QueryResults queries results with pagination and ordering.
func (s *Store) QueryResults(jobID string, orderBy string, offset, limit int) ([]*BlastHit, int, error) {
	// Map order_by to SQL column
	orderCol := "bitscore DESC, evalue ASC"
	switch orderBy {
	case "bitscore":
		orderCol = "bitscore DESC, evalue ASC"
	case "evalue":
		orderCol = "evalue ASC, bitscore DESC"
	case "pident":
		orderCol = "pident DESC, bitscore DESC"
	case "length":
		orderCol = "length DESC, bitscore DESC"
	}

	// Get total count
	var total int
	err := s.db.QueryRow("SELECT COUNT(*) FROM blast_results WHERE job_id = ?", jobID).Scan(&total)
	if err != nil {
		return nil, 0, err
	}

	// Query with pagination
	query := fmt.Sprintf(`
		SELECT dataset_id, gene_id, pident, length, evalue, bitscore
		FROM blast_results
		WHERE job_id = ?
		ORDER BY %s
		LIMIT ? OFFSET ?
	`, orderCol)

	rows, err := s.db.Query(query, jobID, limit, offset)
	if err != nil {
		return nil, 0, err
	}
	defer rows.Close()

	var results []*BlastHit
	for rows.Next() {
		var r BlastHit
		err := rows.Scan(
			&r.DatasetID, &r.GeneID,
			&r.Pident, &r.Length, &r.Evalue, &r.Bitscore,
		)
		if err != nil {
			return nil, 0, err
		}
		results = append(results, &r)
	}

	return results, total, nil
}

// ListQueuedJobs returns all queued jobs (for restart recovery).
func (s *Store) ListQueuedJobs() ([]*BlastJob, error) {
	rows, err := s.db.Query(`
		SELECT job_id, status, params_json, phase, done, total, error, created_at, started_at, finished_at
		FROM blast_jobs WHERE status = ?
		ORDER BY created_at ASC
	`, string(JobStatusQueued))
	if err != nil {
		return nil, err
	}
	defer rows.Close()

	return s.scanJobs(rows)
}

// MarkRunningAsFailed marks all running jobs as failed (for restart recovery).
func (s *Store) MarkRunningAsFailed(errMsg string) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	now := time.Now().Format(time.RFC3339)
	_, err := s.db.Exec(`
		UPDATE blast_jobs SET status = ?, error = ?, finished_at = ?
		WHERE status = ?
	`, string(JobStatusFailed), errMsg, now, string(JobStatusRunning))
	return err
}

// DeleteExpiredJobs deletes jobs older than retentionDays.
func (s *Store) DeleteExpiredJobs(retentionDays int) (int64, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	cutoff := time.Now().AddDate(0, 0, -retentionDays).Format(time.RFC3339)

	// Delete results first (foreign key)
	_, err := s.db.Exec(`
		DELETE FROM blast_results WHERE job_id IN (
			SELECT job_id FROM blast_jobs WHERE finished_at IS NOT NULL AND finished_at < ?
		)
	`, cutoff)
	if err != nil {
		return 0, err
	}

	// Delete jobs
	result, err := s.db.Exec(`
		DELETE FROM blast_jobs WHERE finished_at IS NOT NULL AND finished_at < ?
	`, cutoff)
	if err != nil {
		return 0, err
	}

	return result.RowsAffected()
}

// DeleteJob deletes a job and its results.
func (s *Store) DeleteJob(jobID string) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	// Delete results first
	_, err := s.db.Exec("DELETE FROM blast_results WHERE job_id = ?", jobID)
	if err != nil {
		return err
	}

	_, err = s.db.Exec("DELETE FROM blast_jobs WHERE job_id = ?", jobID)
	return err
}

func (s *Store) scanJobs(rows *sql.Rows) ([]*BlastJob, error) {
	var jobs []*BlastJob
	for rows.Next() {
		var job BlastJob
		var paramsJSON string
		var createdAtStr string
		var startedAtStr, finishedAtStr sql.NullString

		err := rows.Scan(
			&job.ID,
			&job.Status,
			&paramsJSON,
			&job.Progress.Phase,
			&job.Progress.Done,
			&job.Progress.Total,
			&job.Error,
			&createdAtStr,
			&startedAtStr,
			&finishedAtStr,
		)
		if err != nil {
			return nil, err
		}

		if err := json.Unmarshal([]byte(paramsJSON), &job.Params); err != nil {
			return nil, fmt.Errorf("failed to unmarshal params: %w", err)
		}

		job.CreatedAt, _ = time.Parse(time.RFC3339, createdAtStr)
		if startedAtStr.Valid {
			t, _ := time.Parse(time.RFC3339, startedAtStr.String)
			job.StartedAt = &t
		}
		if finishedAtStr.Valid {
			t, _ := time.Parse(time.RFC3339, finishedAtStr.String)
			job.FinishedAt = &t
		}

		jobs = append(jobs, &job)
	}
	return jobs, nil
}

