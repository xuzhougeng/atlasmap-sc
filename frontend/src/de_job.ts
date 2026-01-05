import './styles/main.css';

import { ApiClient, SomaDeJobStatusResponse } from './api/client';
import { ThemeManager } from './components/ThemeManager';

interface DatasetInfo {
    id: string;
    name: string;
}

interface DatasetsResponse {
    default: string;
    datasets: DatasetInfo[];
    title: string;
}

function escapeHtml(text: string): string {
    return text
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;')
        .replace(/'/g, '&#039;');
}

async function fetchDatasets(): Promise<DatasetsResponse> {
    const response = await fetch('/api/datasets');
    if (!response.ok) {
        throw new Error(`Failed to fetch datasets: ${response.statusText}`);
    }
    return response.json();
}

function createDatasetSelector(
    container: HTMLElement,
    datasets: DatasetInfo[],
    currentDataset: string
): void {
    if (datasets.length <= 1) return;

    const selectorContainer = document.createElement('div');
    selectorContainer.className = 'dataset-selector';

    const label = document.createElement('label');
    label.textContent = 'Dataset: ';
    label.htmlFor = 'dataset-select';

    const select = document.createElement('select');
    select.id = 'dataset-select';
    select.className = 'dataset-select';

    for (const ds of datasets) {
        const option = document.createElement('option');
        option.value = ds.id;
        option.textContent = ds.name;
        option.selected = ds.id === currentDataset;
        select.appendChild(option);
    }

    select.addEventListener('change', () => {
        const url = new URL(window.location.href);
        url.searchParams.set('dataset', select.value);
        window.location.href = url.toString();
    });

    selectorContainer.appendChild(label);
    selectorContainer.appendChild(select);
    container.appendChild(selectorContainer);
}

function formatTimestamp(ts?: string): string {
    if (!ts) return '—';
    const d = new Date(ts);
    if (Number.isNaN(d.getTime())) return ts;
    return d.toLocaleString();
}

function renderStatus(root: HTMLElement, info: SomaDeJobStatusResponse, currentUrl: string) {
    const progress = info.progress;
    const hasProgress = progress && (progress.done !== undefined && progress.total !== undefined);
    const progressPercent = hasProgress ? Math.round((progress.done! / progress.total!) * 100) : 0;
    const progressText = hasProgress
        ? `${progress.phase ? progress.phase + ': ' : ''}${progress.done}/${progress.total}`
        : (progress?.phase || '—');

    root.innerHTML = `
        <div class="de-job-status">
            <div class="status-grid">
                <div class="status-section">
                    <div class="status-section-title">Job Identifier</div>
                    <div class="job-id">${escapeHtml(info.job_id)}</div>
                </div>

                <div class="status-section">
                    <div class="status-section-title">Current Status</div>
                    <span class="status-badge ${escapeHtml(info.status)}" id="status-badge">${escapeHtml(info.status)}</span>
                    ${info.status === 'running' || info.status === 'queued' ? `
                        <div class="progress-container">
                            <div class="progress-bar">
                                <div class="progress-fill" style="width: ${progressPercent}%"></div>
                            </div>
                            <div class="progress-text">${escapeHtml(progressText)}</div>
                        </div>
                    ` : ''}
                    ${info.error ? `<div class="error-box">${escapeHtml(info.error)}</div>` : ''}
                </div>

                <div class="status-section">
                    <div class="status-section-title">Timestamps</div>
                    <div class="timestamps">
                        <div class="timestamp-item">
                            <div class="timestamp-label">Created</div>
                            <div class="timestamp-value">${escapeHtml(formatTimestamp(info.created_at))}</div>
                        </div>
                        <div class="timestamp-item">
                            <div class="timestamp-label">Started</div>
                            <div class="timestamp-value">${escapeHtml(formatTimestamp(info.started_at))}</div>
                        </div>
                        <div class="timestamp-item">
                            <div class="timestamp-label">Finished</div>
                            <div class="timestamp-value">${escapeHtml(formatTimestamp(info.finished_at))}</div>
                        </div>
                    </div>
                </div>

                <div class="actions">
                    <button class="btn-action primary" id="btn-refresh" type="button">
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <path d="M23 4v6h-6M1 20v-6h6"/>
                            <path d="M3.51 9a9 9 0 0114.85-3.36L23 10M1 14l4.64 4.36A9 9 0 0020.49 15"/>
                        </svg>
                        Refresh
                    </button>
                    <button class="btn-action" id="btn-copy" type="button">
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <rect x="9" y="9" width="13" height="13" rx="2" ry="2"/>
                            <path d="M5 15H4a2 2 0 01-2-2V4a2 2 0 012-2h9a2 2 0 012 2v1"/>
                        </svg>
                        Copy URL
                    </button>
                </div>

                <div class="url-notice">
                    <div class="url-notice-title">Bookmark this URL</div>
                    <div class="url-notice-url">${escapeHtml(currentUrl)}</div>
                </div>
            </div>
        </div>
    `;
}

async function init() {
    new ThemeManager();

    const root = document.getElementById('de-job-root');
    if (!root) return;

    const url = new URL(window.location.href);
    const dataset = url.searchParams.get('dataset') || '';
    const jobId = url.searchParams.get('job_id') || '';

    if (!dataset || !jobId) {
        root.innerHTML = `<div class="error-box">Missing parameters: requires \`dataset\` and \`job_id\`.</div>`;
        return;
    }

    let datasetsInfo: DatasetsResponse | null = null;
    try {
        datasetsInfo = await fetchDatasets();
    } catch {
        // Optional; page still works without dataset selector.
    }

    if (datasetsInfo?.title) {
        const titleEl = document.querySelector('.header-title');
        if (titleEl) titleEl.textContent = datasetsInfo.title;
        document.title = `${datasetsInfo.title} - DE Job`;
    }

    const datasetInfoContainer = document.getElementById('dataset-info');
    if (datasetInfoContainer && datasetsInfo) {
        createDatasetSelector(datasetInfoContainer, datasetsInfo.datasets, dataset);
    }

    const backBtn = document.getElementById('btn-back') as HTMLAnchorElement | null;
    if (backBtn) {
        backBtn.href = `/?dataset=${encodeURIComponent(dataset)}`;
    }

    const api = new ApiClient(`/d/${encodeURIComponent(dataset)}/api`);

    let pollTimer: number | null = null;
    let inFlight = false;

    const refresh = async () => {
        if (inFlight) return;
        inFlight = true;
        try {
            const info = await api.getSomaDeJob(jobId);
            renderStatus(root, info, window.location.href);

            const refreshBtn = root.querySelector('#btn-refresh') as HTMLButtonElement | null;
            refreshBtn?.addEventListener('click', () => void refresh());

            const copyBtn = root.querySelector('#btn-copy') as HTMLButtonElement | null;
            copyBtn?.addEventListener('click', async () => {
                try {
                    await navigator.clipboard.writeText(window.location.href);
                    if (copyBtn) {
                        copyBtn.innerHTML = `
                            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <polyline points="20 6 9 17 4 12"/>
                            </svg>
                            Copied!
                        `;
                        window.setTimeout(() => {
                            copyBtn.innerHTML = `
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <rect x="9" y="9" width="13" height="13" rx="2" ry="2"/>
                                    <path d="M5 15H4a2 2 0 01-2-2V4a2 2 0 012-2h9a2 2 0 012 2v1"/>
                                </svg>
                                Copy URL
                            `;
                        }, 1500);
                    }
                } catch {
                    if (copyBtn) copyBtn.textContent = 'Copy failed';
                }
            });

            if (info.status === 'completed' || info.status === 'failed' || info.status === 'cancelled') {
                if (pollTimer !== null) {
                    clearInterval(pollTimer);
                    pollTimer = null;
                }
            }
        } catch (e: any) {
            root.innerHTML = `<div class="error-box">Failed to fetch status: ${escapeHtml(e?.message || String(e))}</div>`;
        } finally {
            inFlight = false;
        }
    };

    await refresh();
    pollTimer = window.setInterval(() => void refresh(), 2000);
}

document.addEventListener('DOMContentLoaded', () => {
    void init();
});
