// Package api provides HTTP handlers for the SOMA-Tiles server.
package api

import (
	"context"
	"crypto/rand"
	"encoding/hex"
	"log"
	"sync"
	"time"

	"github.com/soma-tiles/server/internal/destore"
)

// JobManagerConfig contains configuration for the job manager.
type JobManagerConfig struct {
	MaxConcurrent int    // Max concurrent DE jobs (default 1)
	SQLitePath    string // Path to SQLite database
	RetentionDays int    // Days to keep completed jobs (default 7)
	CleanupPeriod time.Duration
}

// JobManager manages DE jobs with SQLite persistence.
type JobManager struct {
	cfg      JobManagerConfig
	store    *destore.Store
	queue    chan string // job IDs
	running  map[string]context.CancelFunc
	mu       sync.Mutex
	wg       sync.WaitGroup
	stopOnce sync.Once
	stopCh   chan struct{}

	// Executor is called to run the actual DE computation.
	Executor func(ctx context.Context, store *destore.Store, jobID string) error
}

// NewJobManager creates a new job manager with SQLite persistence.
func NewJobManager(cfg JobManagerConfig) (*JobManager, error) {
	if cfg.MaxConcurrent <= 0 {
		cfg.MaxConcurrent = 1
	}
	if cfg.RetentionDays <= 0 {
		cfg.RetentionDays = 7
	}
	if cfg.CleanupPeriod <= 0 {
		cfg.CleanupPeriod = 1 * time.Hour
	}

	store, err := destore.NewStore(cfg.SQLitePath)
	if err != nil {
		return nil, err
	}

	jm := &JobManager{
		cfg:     cfg,
		store:   store,
		queue:   make(chan string, 100),
		running: make(map[string]context.CancelFunc),
		stopCh:  make(chan struct{}),
	}
	return jm, nil
}

// Store returns the underlying store for direct access.
func (jm *JobManager) Store() *destore.Store {
	return jm.store
}

// Start starts the worker goroutines and cleanup ticker.
// Also recovers from previous shutdown.
func (jm *JobManager) Start() {
	// Mark any running jobs as failed (server restart)
	if err := jm.store.MarkRunningAsFailed("server restarted"); err != nil {
		log.Printf("[JobManager] failed to mark running jobs as failed: %v", err)
	}

	// Re-queue any queued jobs
	queued, err := jm.store.ListQueuedJobs()
	if err != nil {
		log.Printf("[JobManager] failed to list queued jobs: %v", err)
	} else {
		for _, job := range queued {
			select {
			case jm.queue <- job.ID:
				log.Printf("[JobManager] re-queued job %s", job.ID)
			default:
				log.Printf("[JobManager] queue full, cannot re-queue job %s", job.ID)
			}
		}
	}

	// Start workers
	for i := 0; i < jm.cfg.MaxConcurrent; i++ {
		jm.wg.Add(1)
		go jm.worker()
	}

	// Start cleanup ticker
	go jm.cleaner()
}

// Stop stops all workers gracefully.
func (jm *JobManager) Stop() {
	jm.stopOnce.Do(func() {
		close(jm.stopCh)
		close(jm.queue)
		jm.wg.Wait()
		jm.store.Close()
	})
}

func (jm *JobManager) worker() {
	defer jm.wg.Done()
	for jobID := range jm.queue {
		jm.runJob(jobID)
	}
}

func (jm *JobManager) runJob(jobID string) {
	ctx, cancel := context.WithCancel(context.Background())

	jm.mu.Lock()
	jm.running[jobID] = cancel
	jm.mu.Unlock()

	defer func() {
		jm.mu.Lock()
		delete(jm.running, jobID)
		jm.mu.Unlock()
	}()

	// Mark as running
	if err := jm.store.UpdateJobStarted(jobID); err != nil {
		log.Printf("[JobManager] failed to update job %s as started: %v", jobID, err)
		return
	}

	var execErr error
	if jm.Executor != nil {
		execErr = jm.Executor(ctx, jm.store, jobID)
	}

	// Update final status
	if ctx.Err() == context.Canceled {
		jm.store.UpdateJobStatus(jobID, destore.JobStatusCancelled, "cancelled by user")
	} else if execErr != nil {
		jm.store.UpdateJobStatus(jobID, destore.JobStatusFailed, execErr.Error())
	} else {
		jm.store.UpdateJobStatus(jobID, destore.JobStatusCompleted, "")
	}
}

func (jm *JobManager) cleaner() {
	ticker := time.NewTicker(jm.cfg.CleanupPeriod)
	defer ticker.Stop()
	for {
		select {
		case <-jm.stopCh:
			return
		case <-ticker.C:
			jm.cleanup()
		}
	}
}

func (jm *JobManager) cleanup() {
	deleted, err := jm.store.DeleteExpiredJobs(jm.cfg.RetentionDays)
	if err != nil {
		log.Printf("[JobManager] cleanup error: %v", err)
	} else if deleted > 0 {
		log.Printf("[JobManager] cleaned up %d expired jobs", deleted)
	}
}

// Submit creates a new job and enqueues it for execution.
func (jm *JobManager) Submit(params destore.DEJobParams) (*destore.DEJob, error) {
	id := generateJobID()
	job := &destore.DEJob{
		ID:        id,
		DatasetID: params.DatasetID,
		Status:    destore.JobStatusQueued,
		Params:    params,
		CreatedAt: time.Now(),
	}

	if err := jm.store.CreateJob(job); err != nil {
		return nil, err
	}

	select {
	case jm.queue <- id:
	default:
		// Queue full; mark as failed immediately
		jm.store.UpdateJobStatus(id, destore.JobStatusFailed, "job queue is full; try again later")
	}

	return job, nil
}

// Get returns a job by ID.
func (jm *JobManager) Get(id string) *destore.DEJob {
	job, err := jm.store.GetJob(id)
	if err != nil {
		log.Printf("[JobManager] error getting job %s: %v", id, err)
		return nil
	}
	return job
}

// Cancel attempts to cancel a running job.
func (jm *JobManager) Cancel(id string) bool {
	jm.mu.Lock()
	cancel, ok := jm.running[id]
	jm.mu.Unlock()

	if ok && cancel != nil {
		cancel()
		return true
	}

	// If not running, try to mark as cancelled in DB
	job, err := jm.store.GetJob(id)
	if err != nil || job == nil {
		return false
	}
	if job.Status == destore.JobStatusQueued {
		jm.store.UpdateJobStatus(id, destore.JobStatusCancelled, "cancelled before start")
		return true
	}
	return false
}

// Delete deletes a job and its results.
func (jm *JobManager) Delete(id string) error {
	return jm.store.DeleteJob(id)
}

func generateJobID() string {
	b := make([]byte, 8)
	rand.Read(b)
	return hex.EncodeToString(b)
}
