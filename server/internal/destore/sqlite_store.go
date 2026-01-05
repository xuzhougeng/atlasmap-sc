// Package destore provides persistent storage for DE job state and results using SQLite.
package destore

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

// JobStatus represents the current state of a DE job.
type JobStatus string

const (
	JobStatusQueued    JobStatus = "queued"
	JobStatusRunning   JobStatus = "running"
	JobStatusCompleted JobStatus = "completed"
	JobStatusFailed    JobStatus = "failed"
	JobStatusCancelled JobStatus = "cancelled"
)

// DEJobParams contains the parameters for a DE job.
type DEJobParams struct {
	DatasetID        string   `json:"dataset_id"`
	Groupby          string   `json:"groupby"`
	Group1           []string `json:"group1"`
	Group2           []string `json:"group2"`
	Tests            []string `json:"tests"`
	MaxCellsPerGroup int      `json:"max_cells_per_group"`
	Seed             int      `json:"seed"`
	Limit            int      `json:"limit"`
}

// DEJobProgress represents the progress of a DE job.
type DEJobProgress struct {
	Phase string `json:"phase"`
	Done  int    `json:"done"`
	Total int    `json:"total"`
}

// DEJob represents a differential expression job.
type DEJob struct {
	ID         string        `json:"job_id"`
	DatasetID  string        `json:"dataset_id"`
	Status     JobStatus     `json:"status"`
	Params     DEJobParams   `json:"params"`
	Progress   DEJobProgress `json:"progress"`
	CreatedAt  time.Time     `json:"created_at"`
	StartedAt  *time.Time    `json:"started_at,omitempty"`
	FinishedAt *time.Time    `json:"finished_at,omitempty"`
	N1         int           `json:"n1"`
	N2         int           `json:"n2"`
	Error      string        `json:"error,omitempty"`
}

// DEGeneResult contains the DE result for a single gene.
type DEGeneResult struct {
	Gene       string  `json:"gene"`
	GeneJoinID int64   `json:"gene_joinid"`
	Mean1      float64 `json:"mean1"`
	Mean2      float64 `json:"mean2"`
	Pct1       float64 `json:"pct1"`
	Pct2       float64 `json:"pct2"`
	Log2FC     float64 `json:"log2fc"`
	PTtest     float64 `json:"p_ttest"`
	FDRTtest   float64 `json:"fdr_ttest"`
	PRanksum   float64 `json:"p_ranksum"`
	FDRRanksum float64 `json:"fdr_ranksum"`
}

// Store provides persistent storage for DE jobs using SQLite.
type Store struct {
	db *sql.DB
	mu sync.Mutex
}

// NewStore creates a new SQLite-based DE store.
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
	CREATE TABLE IF NOT EXISTS de_jobs (
		job_id TEXT PRIMARY KEY,
		dataset_id TEXT NOT NULL,
		status TEXT NOT NULL,
		params_json TEXT NOT NULL,
		phase TEXT DEFAULT '',
		done INTEGER DEFAULT 0,
		total INTEGER DEFAULT 0,
		n1 INTEGER DEFAULT 0,
		n2 INTEGER DEFAULT 0,
		error TEXT DEFAULT '',
		created_at TEXT NOT NULL,
		started_at TEXT,
		finished_at TEXT
	);

	CREATE INDEX IF NOT EXISTS idx_de_jobs_dataset ON de_jobs(dataset_id);
	CREATE INDEX IF NOT EXISTS idx_de_jobs_status ON de_jobs(status);
	CREATE INDEX IF NOT EXISTS idx_de_jobs_finished ON de_jobs(finished_at);

	CREATE TABLE IF NOT EXISTS de_results (
		id INTEGER PRIMARY KEY AUTOINCREMENT,
		job_id TEXT NOT NULL,
		gene TEXT NOT NULL,
		gene_joinid INTEGER NOT NULL,
		mean1 REAL NOT NULL,
		mean2 REAL NOT NULL,
		pct1 REAL NOT NULL,
		pct2 REAL NOT NULL,
		log2fc REAL NOT NULL,
		p_ttest REAL NOT NULL,
		fdr_ttest REAL NOT NULL,
		p_ranksum REAL NOT NULL,
		fdr_ranksum REAL NOT NULL,
		FOREIGN KEY (job_id) REFERENCES de_jobs(job_id) ON DELETE CASCADE
	);

	CREATE INDEX IF NOT EXISTS idx_de_results_job ON de_results(job_id);
	CREATE INDEX IF NOT EXISTS idx_de_results_job_fdr_ranksum ON de_results(job_id, fdr_ranksum);
	CREATE INDEX IF NOT EXISTS idx_de_results_job_fdr_ttest ON de_results(job_id, fdr_ttest);
	`
	_, err := s.db.Exec(schema)
	return err
}

// CreateJob creates a new job record with status=queued.
func (s *Store) CreateJob(job *DEJob) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	paramsJSON, err := json.Marshal(job.Params)
	if err != nil {
		return fmt.Errorf("failed to marshal params: %w", err)
	}

	_, err = s.db.Exec(`
		INSERT INTO de_jobs (job_id, dataset_id, status, params_json, phase, done, total, n1, n2, error, created_at, started_at, finished_at)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	`,
		job.ID,
		job.Params.DatasetID,
		string(job.Status),
		string(paramsJSON),
		job.Progress.Phase,
		job.Progress.Done,
		job.Progress.Total,
		job.N1,
		job.N2,
		job.Error,
		job.CreatedAt.Format(time.RFC3339),
		nil,
		nil,
	)
	return err
}

// GetJob retrieves a job by ID.
func (s *Store) GetJob(jobID string) (*DEJob, error) {
	row := s.db.QueryRow(`
		SELECT job_id, dataset_id, status, params_json, phase, done, total, n1, n2, error, created_at, started_at, finished_at
		FROM de_jobs WHERE job_id = ?
	`, jobID)

	var job DEJob
	var paramsJSON string
	var createdAtStr string
	var startedAtStr, finishedAtStr sql.NullString

	err := row.Scan(
		&job.ID,
		&job.DatasetID,
		&job.Status,
		&paramsJSON,
		&job.Progress.Phase,
		&job.Progress.Done,
		&job.Progress.Total,
		&job.N1,
		&job.N2,
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
		UPDATE de_jobs SET status = ?, error = ?, finished_at = COALESCE(?, finished_at)
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
		UPDATE de_jobs SET status = ?, started_at = ?
		WHERE job_id = ?
	`, string(JobStatusRunning), now, jobID)
	return err
}

// UpdateJobProgress updates the progress fields.
func (s *Store) UpdateJobProgress(jobID string, phase string, done, total int) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	_, err := s.db.Exec(`
		UPDATE de_jobs SET phase = ?, done = ?, total = ?
		WHERE job_id = ?
	`, phase, done, total, jobID)
	return err
}

// UpdateJobCounts updates n1 and n2.
func (s *Store) UpdateJobCounts(jobID string, n1, n2 int) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	_, err := s.db.Exec(`
		UPDATE de_jobs SET n1 = ?, n2 = ?
		WHERE job_id = ?
	`, n1, n2, jobID)
	return err
}

// InsertResults inserts gene results in a batch transaction.
func (s *Store) InsertResults(jobID string, results []*DEGeneResult) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	tx, err := s.db.Begin()
	if err != nil {
		return err
	}
	defer tx.Rollback()

	stmt, err := tx.Prepare(`
		INSERT INTO de_results (job_id, gene, gene_joinid, mean1, mean2, pct1, pct2, log2fc, p_ttest, fdr_ttest, p_ranksum, fdr_ranksum)
		VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	`)
	if err != nil {
		return err
	}
	defer stmt.Close()

	for _, r := range results {
		_, err := stmt.Exec(
			jobID, r.Gene, r.GeneJoinID,
			r.Mean1, r.Mean2, r.Pct1, r.Pct2, r.Log2FC,
			r.PTtest, r.FDRTtest, r.PRanksum, r.FDRRanksum,
		)
		if err != nil {
			return err
		}
	}

	return tx.Commit()
}

// QueryResults queries results with pagination and ordering.
func (s *Store) QueryResults(jobID string, orderBy string, offset, limit int) ([]*DEGeneResult, int, error) {
	// Map order_by to SQL column
	orderCol := "fdr_ranksum ASC, ABS(log2fc) DESC"
	switch orderBy {
	case "fdr_ranksum":
		orderCol = "fdr_ranksum ASC, ABS(log2fc) DESC"
	case "fdr_ttest":
		orderCol = "fdr_ttest ASC, ABS(log2fc) DESC"
	case "p_ranksum":
		orderCol = "p_ranksum ASC, ABS(log2fc) DESC"
	case "p_ttest":
		orderCol = "p_ttest ASC, ABS(log2fc) DESC"
	case "abs_log2fc":
		orderCol = "ABS(log2fc) DESC, fdr_ranksum ASC"
	}

	// Get total count
	var total int
	err := s.db.QueryRow("SELECT COUNT(*) FROM de_results WHERE job_id = ?", jobID).Scan(&total)
	if err != nil {
		return nil, 0, err
	}

	// Query with pagination
	query := fmt.Sprintf(`
		SELECT gene, gene_joinid, mean1, mean2, pct1, pct2, log2fc, p_ttest, fdr_ttest, p_ranksum, fdr_ranksum
		FROM de_results
		WHERE job_id = ?
		ORDER BY %s
		LIMIT ? OFFSET ?
	`, orderCol)

	rows, err := s.db.Query(query, jobID, limit, offset)
	if err != nil {
		return nil, 0, err
	}
	defer rows.Close()

	var results []*DEGeneResult
	for rows.Next() {
		var r DEGeneResult
		err := rows.Scan(
			&r.Gene, &r.GeneJoinID,
			&r.Mean1, &r.Mean2, &r.Pct1, &r.Pct2, &r.Log2FC,
			&r.PTtest, &r.FDRTtest, &r.PRanksum, &r.FDRRanksum,
		)
		if err != nil {
			return nil, 0, err
		}
		results = append(results, &r)
	}

	return results, total, nil
}

// ListJobsByDataset returns all jobs for a dataset.
func (s *Store) ListJobsByDataset(datasetID string) ([]*DEJob, error) {
	rows, err := s.db.Query(`
		SELECT job_id, dataset_id, status, params_json, phase, done, total, n1, n2, error, created_at, started_at, finished_at
		FROM de_jobs WHERE dataset_id = ?
		ORDER BY created_at DESC
	`, datasetID)
	if err != nil {
		return nil, err
	}
	defer rows.Close()

	return s.scanJobs(rows)
}

// ListQueuedJobs returns all queued jobs (for restart recovery).
func (s *Store) ListQueuedJobs() ([]*DEJob, error) {
	rows, err := s.db.Query(`
		SELECT job_id, dataset_id, status, params_json, phase, done, total, n1, n2, error, created_at, started_at, finished_at
		FROM de_jobs WHERE status = ?
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
		UPDATE de_jobs SET status = ?, error = ?, finished_at = ?
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
		DELETE FROM de_results WHERE job_id IN (
			SELECT job_id FROM de_jobs WHERE finished_at IS NOT NULL AND finished_at < ?
		)
	`, cutoff)
	if err != nil {
		return 0, err
	}

	// Delete jobs
	result, err := s.db.Exec(`
		DELETE FROM de_jobs WHERE finished_at IS NOT NULL AND finished_at < ?
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
	_, err := s.db.Exec("DELETE FROM de_results WHERE job_id = ?", jobID)
	if err != nil {
		return err
	}

	_, err = s.db.Exec("DELETE FROM de_jobs WHERE job_id = ?", jobID)
	return err
}

func (s *Store) scanJobs(rows *sql.Rows) ([]*DEJob, error) {
	var jobs []*DEJob
	for rows.Next() {
		var job DEJob
		var paramsJSON string
		var createdAtStr string
		var startedAtStr, finishedAtStr sql.NullString

		err := rows.Scan(
			&job.ID,
			&job.DatasetID,
			&job.Status,
			&paramsJSON,
			&job.Progress.Phase,
			&job.Progress.Done,
			&job.Progress.Total,
			&job.N1,
			&job.N2,
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

