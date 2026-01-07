// BLASTP Result Page Entry Point

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

interface BlastHit {
    dataset_id: string;
    gene_id: string;
    pident: number;
    length: number;
    evalue: number;
    bitscore: number;
}

interface BlastJobResultResponse {
    params: {
        sequence: string;
        max_hits: number;
        evalue: number;
        datasets?: string[];
    };
    total: number;
    offset: number;
    limit: number;
    order_by: string;
    items: BlastHit[];
}

async function fetchDatasets(): Promise<DatasetsResponse> {
    const response = await fetch('/api/datasets');
    if (!response.ok) {
        throw new Error(`Failed to fetch datasets: ${response.statusText}`);
    }
    return response.json();
}

function escapeHtml(text: string): string {
    return text
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
}

function formatNumber(v: number, digits: number = 4): string {
    if (!Number.isFinite(v)) return '-';
    const abs = Math.abs(v);
    if (abs !== 0 && abs < 1e-4) return v.toExponential(2);
    if (abs >= 1e6) return v.toExponential(2);
    return v.toFixed(digits);
}

function setBadgeClass(el: HTMLElement, status: BlastJobStatus): void {
    el.classList.remove('queued', 'running', 'completed', 'failed', 'cancelled');
    el.classList.add(status);
}

function renderTableBody(items: BlastHit[]): string {
    if (!items.length) {
        return `<tr><td colspan="6" style="color:var(--text-secondary); padding:16px; text-align:center;">No hits found</td></tr>`;
    }

    return items
        .map((hit) => {
            const geneLink = `/?dataset=${encodeURIComponent(hit.dataset_id)}&gene=${encodeURIComponent(hit.gene_id)}`;
            return `
                <tr>
                    <td><span class="blastp-dataset-badge">${escapeHtml(hit.dataset_id)}</span></td>
                    <td><a href="${geneLink}" class="gene-link" title="View expression pattern">${escapeHtml(hit.gene_id)}</a></td>
                    <td>${formatNumber(hit.pident, 2)}%</td>
                    <td>${hit.length}</td>
                    <td>${formatNumber(hit.evalue, 4)}</td>
                    <td>${formatNumber(hit.bitscore, 1)}</td>
                </tr>
            `;
        })
        .join('');
}

function renderResult(root: HTMLElement, result: BlastJobResultResponse): void {
    const params = result.params;
    const seqPreview = params.sequence.length > 50 
        ? params.sequence.substring(0, 50) + '...' 
        : params.sequence;

    // Summary
    root.querySelector('#summary')!.innerHTML = `
        <div class="de-result-summary-item">
            <span class="label">Query</span>
            <span class="value mono" title="${escapeHtml(params.sequence)}">${escapeHtml(seqPreview)}</span>
        </div>
        <div class="de-result-summary-item">
            <span class="label">E-value Cutoff</span>
            <span class="value">${formatNumber(params.evalue, 4)}</span>
        </div>
        <div class="de-result-summary-item">
            <span class="label">Max Hits</span>
            <span class="value">${params.max_hits}</span>
        </div>
        <div class="de-result-summary-divider"></div>
        <div class="de-result-summary-item">
            <span class="label">Total Hits</span>
            <span class="value">${result.total.toLocaleString()}</span>
        </div>
    `;

    // Table body
    const tableBody = root.querySelector('#results-body')!;
    tableBody.innerHTML = renderTableBody(result.items);

    // Pagination info
    const pageInfoEl = root.querySelector('#page-info') as HTMLSpanElement;
    if (pageInfoEl) {
        const start = result.total === 0 ? 0 : result.offset + 1;
        const end = Math.min(result.total, result.offset + result.limit);
        pageInfoEl.textContent = `${start}-${end} / ${result.total}`;
    }

    // Order by
    const orderByEl = root.querySelector('#order-by') as HTMLSelectElement;
    if (orderByEl) orderByEl.value = result.order_by;

    // Pagination buttons
    const prevBtn = root.querySelector('#btn-prev') as HTMLButtonElement;
    const nextBtn = root.querySelector('#btn-next') as HTMLButtonElement;
    prevBtn.disabled = result.offset <= 0;
    nextBtn.disabled = result.offset + result.limit >= result.total;
}

function renderBase(root: HTMLElement, jobId: string): void {
    root.innerHTML = `
        <!-- Status Strip -->
        <div class="de-result-status-strip">
            <div class="de-result-status-left">
                <div class="de-result-job-id">
                    <span class="label">Job ID</span>
                    <span class="value">${escapeHtml(jobId)}</span>
                </div>
                <span class="de-result-status-badge queued" id="status-badge">loading</span>
            </div>
            <div class="de-result-status-right">
                <button class="de-result-btn" id="btn-refresh" type="button">
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <path d="M23 4v6h-6M1 20v-6h6"/>
                        <path d="M3.51 9a9 9 0 0 1 14.85-3.36L23 10M1 14l4.64 4.36A9 9 0 0 0 20.49 15"/>
                    </svg>
                    Refresh
                </button>
                <button class="de-result-btn" id="btn-copy" type="button">
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
                        <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>
                    </svg>
                    Copy Link
                </button>
            </div>
        </div>

        <!-- Summary Strip -->
        <div class="de-result-summary-strip" id="summary"></div>

        <!-- Error Container -->
        <div id="error" class="de-result-error" style="display:none;"></div>

        <!-- Table Panel -->
        <div class="de-result-panel de-result-table-panel blastp-table-panel">
            <div class="de-result-panel-header">
                <span class="de-result-panel-title">BLAST Hits</span>
                <div class="de-result-panel-controls">
                    <div class="de-result-control-group">
                        <span class="label">Sort By</span>
                        <select id="order-by" class="de-result-select">
                            <option value="bitscore">bitscore</option>
                            <option value="evalue">evalue</option>
                            <option value="pident">pident</option>
                            <option value="length">length</option>
                        </select>
                    </div>
                    <div class="de-result-control-group">
                        <span class="label">Limit</span>
                        <input id="limit" class="de-result-input" type="number" min="1" max="500" value="50" />
                    </div>
                    <button class="de-result-btn-sm" id="btn-apply" type="button">Apply</button>
                </div>
            </div>
            <div class="de-result-table-body">
                <table class="de-result-table blastp-result-table">
                    <thead>
                        <tr>
                            <th>Dataset</th>
                            <th>Gene</th>
                            <th>Identity</th>
                            <th>Length</th>
                            <th>E-value</th>
                            <th>Bitscore</th>
                        </tr>
                    </thead>
                    <tbody id="results-body"></tbody>
                </table>
            </div>
            <div class="de-result-table-footer">
                <div class="de-result-pagination">
                    <button class="de-result-btn-sm" id="btn-prev" type="button">Prev</button>
                    <span id="page-info" class="de-result-pagination-info"></span>
                    <button class="de-result-btn-sm" id="btn-next" type="button">Next</button>
                </div>
            </div>
        </div>
    `;

    const copyBtn = root.querySelector('#btn-copy') as HTMLButtonElement | null;
    copyBtn?.addEventListener('click', async () => {
        try {
            await navigator.clipboard.writeText(window.location.href);
            copyBtn.textContent = 'Copied!';
            setTimeout(() => (copyBtn.innerHTML = `
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
                    <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>
                </svg>
                Copy Link
            `), 1500);
        } catch {
            copyBtn.textContent = 'Failed';
            setTimeout(() => (copyBtn.textContent = 'Copy Link'), 1500);
        }
    });
}

async function init(): Promise<void> {
    new ThemeManager();

    const root = document.getElementById('blastp-result-root');
    if (!root) return;

    const url = new URL(window.location.href);
    const jobId = url.searchParams.get('job_id') || '';

    if (!jobId) {
        root.innerHTML = `<div class="error-box">Missing job_id parameter</div>`;
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
        document.title = `${datasetsInfo.title} - BLASTP Result`;
    }

    // Update back button and job link
    const backBtn = document.getElementById('btn-back') as HTMLAnchorElement | null;
    if (backBtn) {
        backBtn.href = '/';
    }

    const jobBtn = document.getElementById('btn-job') as HTMLAnchorElement | null;
    if (jobBtn) {
        jobBtn.href = `/blastp_job.html?job_id=${encodeURIComponent(jobId)}`;
    }

    renderBase(root, jobId);

    const badge = root.querySelector('#status-badge') as HTMLElement;
    const refreshBtn = root.querySelector('#btn-refresh') as HTMLButtonElement;
    const orderByEl = root.querySelector('#order-by') as HTMLSelectElement;
    const limitEl = root.querySelector('#limit') as HTMLInputElement;
    const applyBtn = root.querySelector('#btn-apply') as HTMLButtonElement;
    const prevBtn = root.querySelector('#btn-prev') as HTMLButtonElement;
    const nextBtn = root.querySelector('#btn-next') as HTMLButtonElement;
    const errorEl = root.querySelector('#error') as HTMLDivElement;

    let lastResult: BlastJobResultResponse | null = null;
    let inFlight = false;

    const renderError = (message: string | null) => {
        if (!message) {
            errorEl.innerHTML = '';
            errorEl.style.display = 'none';
            return;
        }
        errorEl.style.display = 'block';
        errorEl.innerHTML = `<div class="de-result-error-box">${escapeHtml(message)}</div>`;
    };

    const getQueryState = () => {
        const u = new URL(window.location.href);
        const offset = Math.max(0, parseInt(u.searchParams.get('offset') || '0', 10) || 0);
        const limit = Math.max(1, Math.min(500, parseInt(u.searchParams.get('limit') || limitEl.value || '50', 10) || 50));
        const order_by = u.searchParams.get('order_by') || orderByEl.value || 'bitscore';
        return { offset, limit, order_by };
    };

    const setQueryState = (next: { offset?: number; limit?: number; order_by?: string }) => {
        const u = new URL(window.location.href);
        if (next.offset !== undefined) u.searchParams.set('offset', String(next.offset));
        if (next.limit !== undefined) u.searchParams.set('limit', String(next.limit));
        if (next.order_by !== undefined) u.searchParams.set('order_by', next.order_by);
        window.history.replaceState({}, '', u.toString());
    };

    const refreshStatus = async (): Promise<BlastJobStatus> => {
        const response = await fetch(`/api/blastp/jobs/${encodeURIComponent(jobId)}`);
        if (!response.ok) {
            throw new Error(await response.text() || response.statusText);
        }
        const info: BlastJobStatusResponse = await response.json();
        badge.textContent = info.status;
        setBadgeClass(badge, info.status);
        return info.status;
    };

    const loadResults = async () => {
        const { offset, limit, order_by } = getQueryState();
        setQueryState({ offset, limit, order_by });
        limitEl.value = String(limit);
        orderByEl.value = order_by;

        const status = await refreshStatus();
        if (status !== 'completed') {
            renderError(`Job is not completed yet (status: ${status}). Please check the Job Status page.`);
            return;
        }

        const response = await fetch(`/api/blastp/jobs/${encodeURIComponent(jobId)}/result?offset=${offset}&limit=${limit}&order_by=${order_by}`);
        if (!response.ok) {
            throw new Error(await response.text() || response.statusText);
        }
        const result: BlastJobResultResponse = await response.json();
        lastResult = result;
        renderError(null);
        renderResult(root, result);
    };

    const refreshAll = async () => {
        if (inFlight) return;
        inFlight = true;
        refreshBtn.disabled = true;
        try {
            await loadResults();
        } catch (e: any) {
            console.error(e);
            renderError(e?.message || String(e));
        } finally {
            refreshBtn.disabled = false;
            inFlight = false;
        }
    };

    refreshBtn.addEventListener('click', () => void refreshAll());

    applyBtn.addEventListener('click', () => {
        const limit = Math.max(1, Math.min(500, parseInt(limitEl.value || '50', 10) || 50));
        const order_by = orderByEl.value;
        setQueryState({ offset: 0, limit, order_by });
        void refreshAll();
    });

    orderByEl.addEventListener('change', () => {
        const limit = Math.max(1, Math.min(500, parseInt(limitEl.value || '50', 10) || 50));
        setQueryState({ offset: 0, limit, order_by: orderByEl.value });
        void refreshAll();
    });

    prevBtn.addEventListener('click', () => {
        if (!lastResult) return;
        const nextOffset = Math.max(0, lastResult.offset - lastResult.limit);
        setQueryState({ offset: nextOffset });
        void refreshAll();
    });

    nextBtn.addEventListener('click', () => {
        if (!lastResult) return;
        const nextOffset = Math.min(lastResult.total, lastResult.offset + lastResult.limit);
        setQueryState({ offset: nextOffset });
        void refreshAll();
    });

    const state = getQueryState();
    limitEl.value = String(state.limit);
    orderByEl.value = state.order_by;

    await refreshAll();
}

document.addEventListener('DOMContentLoaded', () => {
    void init();
});

