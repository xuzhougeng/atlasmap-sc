// BLASTP Job Status Page Entry Point

import './styles/main.css';

import { ThemeManager } from './components/ThemeManager';

interface DatasetsResponse {
    default: string;
    datasets: { id: string; name: string }[];
    title: string;
}

type BlastJobStatus = 'queued' | 'running' | 'completed' | 'failed' | 'cancelled';

interface BlastJobStatusResponse {
    job_id: string;
    status: BlastJobStatus;
    created_at: string;
    started_at?: string;
    finished_at?: string;
    progress?: {
        phase: string;
        done: number;
        total: number;
    };
    error?: string;
}

async function fetchDatasets(): Promise<DatasetsResponse> {
    const response = await fetch('/api/datasets');
    if (!response.ok) {
        throw new Error(`Failed to fetch datasets: ${response.statusText}`);
    }
    return response.json();
}

function showError(message: string): void {
    const root = document.getElementById('blastp-job-root');
    if (root) {
        root.innerHTML = `<div class="error-box">${message.replace(/</g, '&lt;')}</div>`;
    }
}

function escapeHtml(text: string): string {
    return text
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;');
}

function formatDate(dateStr: string | undefined): string {
    if (!dateStr) return '-';
    try {
        return new Date(dateStr).toLocaleString();
    } catch {
        return dateStr;
    }
}

function renderStatus(root: HTMLElement, jobId: string, status: BlastJobStatusResponse): void {
    const progressHtml = status.progress && status.progress.total > 0
        ? `<div class="de-job-progress">
             <div class="de-job-progress-bar" style="width: ${(status.progress.done / status.progress.total * 100).toFixed(1)}%"></div>
           </div>
           <div class="de-job-progress-text">${status.progress.phase}: ${status.progress.done}/${status.progress.total}</div>`
        : '';

    root.innerHTML = `
        <div class="de-job-status">
            <div class="de-job-row">
                <span class="de-job-label">Job ID</span>
                <span class="de-job-value mono">${escapeHtml(jobId)}</span>
            </div>
            <div class="de-job-row">
                <span class="de-job-label">Status</span>
                <span class="de-job-badge ${status.status}" id="status-badge">${status.status}</span>
            </div>
            ${progressHtml}
            <div class="de-job-row">
                <span class="de-job-label">Created</span>
                <span class="de-job-value">${formatDate(status.created_at)}</span>
            </div>
            <div class="de-job-row">
                <span class="de-job-label">Started</span>
                <span class="de-job-value">${formatDate(status.started_at)}</span>
            </div>
            <div class="de-job-row">
                <span class="de-job-label">Finished</span>
                <span class="de-job-value">${formatDate(status.finished_at)}</span>
            </div>
            ${status.error ? `
            <div class="de-job-row error">
                <span class="de-job-label">Error</span>
                <span class="de-job-value error">${escapeHtml(status.error)}</span>
            </div>
            ` : ''}
        </div>
        
        <div class="de-job-actions">
            <button class="btn-primary" id="btn-refresh">Refresh Status</button>
            ${status.status === 'completed' ? `
            <a class="btn-primary" id="btn-results" href="/blastp_result.html?job_id=${encodeURIComponent(jobId)}">View Results</a>
            ` : ''}
            ${status.status === 'queued' || status.status === 'running' ? `
            <button class="btn-secondary" id="btn-cancel">Cancel Job</button>
            ` : ''}
            <button class="btn-secondary" id="btn-copy">Copy Link</button>
        </div>
    `;

    // Refresh button
    const refreshBtn = root.querySelector('#btn-refresh');
    refreshBtn?.addEventListener('click', () => void refreshStatus(root, jobId));

    // Cancel button
    const cancelBtn = root.querySelector('#btn-cancel');
    cancelBtn?.addEventListener('click', () => void cancelJob(root, jobId));

    // Copy button
    const copyBtn = root.querySelector('#btn-copy') as HTMLButtonElement | null;
    copyBtn?.addEventListener('click', async () => {
        try {
            await navigator.clipboard.writeText(window.location.href);
            copyBtn.textContent = 'Copied!';
            setTimeout(() => (copyBtn.textContent = 'Copy Link'), 1500);
        } catch {
            copyBtn.textContent = 'Failed';
            setTimeout(() => (copyBtn.textContent = 'Copy Link'), 1500);
        }
    });
}

async function refreshStatus(root: HTMLElement, jobId: string): Promise<void> {
    try {
        const response = await fetch(`/api/blastp/jobs/${encodeURIComponent(jobId)}`);
        if (!response.ok) {
            throw new Error(await response.text() || response.statusText);
        }
        const status: BlastJobStatusResponse = await response.json();
        renderStatus(root, jobId, status);

        // Auto-refresh if still running
        if (status.status === 'queued' || status.status === 'running') {
            setTimeout(() => void refreshStatus(root, jobId), 3000);
        }
    } catch (e: any) {
        console.error(e);
        showError('Failed to fetch job status: ' + (e?.message || e));
    }
}

async function cancelJob(root: HTMLElement, jobId: string): Promise<void> {
    try {
        const response = await fetch(`/api/blastp/jobs/${encodeURIComponent(jobId)}`, {
            method: 'DELETE',
        });
        if (!response.ok) {
            throw new Error(await response.text() || response.statusText);
        }
        // Refresh to show cancelled status
        await refreshStatus(root, jobId);
    } catch (e: any) {
        console.error(e);
        alert('Failed to cancel job: ' + (e?.message || e));
    }
}

async function init(): Promise<void> {
    new ThemeManager();

    const root = document.getElementById('blastp-job-root');
    if (!root) return;

    const url = new URL(window.location.href);
    const jobId = url.searchParams.get('job_id') || '';

    if (!jobId) {
        showError('Missing job_id parameter');
        return;
    }

    // Load datasets for title
    let datasetsInfo: DatasetsResponse | null = null;
    try {
        datasetsInfo = await fetchDatasets();
    } catch {
        // Optional
    }

    if (datasetsInfo?.title) {
        const titleEl = document.querySelector('.header-title');
        if (titleEl) titleEl.textContent = datasetsInfo.title;
        document.title = `${datasetsInfo.title} - BLASTP Job`;
    }

    // Load initial status
    await refreshStatus(root, jobId);
}

document.addEventListener('DOMContentLoaded', () => {
    void init();
});

