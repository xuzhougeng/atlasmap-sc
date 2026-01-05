import './styles/main.css';

import {
    ApiClient,
    SomaDeJobResultResponse,
    SomaDeJobStatus,
    SomaDeJobStatusResponse,
    SomaDeResultItem,
    SomaDeResultOrderBy,
} from './api/client';
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

function setBadgeClass(el: HTMLElement, status: SomaDeJobStatus) {
    el.classList.remove('queued', 'running', 'completed', 'failed', 'cancelled');
    el.classList.add(status);
}

function formatNumber(v: number, digits: number = 4): string {
    if (!Number.isFinite(v)) return '-';
    const abs = Math.abs(v);
    if (abs !== 0 && abs < 1e-4) return v.toExponential(2);
    if (abs >= 1e6) return v.toExponential(2);
    return v.toFixed(digits);
}

function formatPct(v: number): string {
    if (!Number.isFinite(v)) return '-';
    return `${(v * 100).toFixed(2)}%`;
}

function renderTableBody(items: SomaDeResultItem[]): string {
    if (!items.length) {
        return `<tr><td colspan="10" style="color:var(--text-secondary); padding:12px;">No results</td></tr>`;
    }

    return items
        .map((it) => {
            return `
                <tr>
                    <td class="mono">${escapeHtml(it.gene)}</td>
                    <td>${formatNumber(it.log2fc, 3)}</td>
                    <td>${formatNumber(it.mean1, 4)}</td>
                    <td>${formatNumber(it.mean2, 4)}</td>
                    <td>${formatPct(it.pct1)}</td>
                    <td>${formatPct(it.pct2)}</td>
                    <td>${formatNumber(it.p_ttest, 4)}</td>
                    <td>${formatNumber(it.fdr_ttest, 4)}</td>
                    <td>${formatNumber(it.p_ranksum, 4)}</td>
                    <td>${formatNumber(it.fdr_ranksum, 4)}</td>
                </tr>
            `;
        })
        .join('');
}

function renderResult(root: HTMLElement, result: SomaDeJobResultResponse, statusUrl: string) {
    const params = result.params;
    const group2Label = params.group2?.length ? params.group2.join(', ') : '(one-vs-rest)';
    root.querySelector('#summary')!.innerHTML = `
        <div><span style="color:var(--text-secondary);">groupby:</span> <span class="mono">${escapeHtml(params.groupby)}</span></div>
        <div><span style="color:var(--text-secondary);">group1:</span> <span class="mono">${escapeHtml(params.group1.join(', '))}</span></div>
        <div><span style="color:var(--text-secondary);">group2:</span> <span class="mono">${escapeHtml(group2Label)}</span></div>
        <div><span style="color:var(--text-secondary);">tests:</span> <span class="mono">${escapeHtml((params.tests || []).join(', '))}</span></div>
        <div><span style="color:var(--text-secondary);">n1/n2:</span> <span class="mono">${result.n1} / ${result.n2}</span></div>
        <div><span style="color:var(--text-secondary);">total genes:</span> <span class="mono">${result.total}</span></div>
        <div style="margin-top:8px; color:var(--text-secondary);">任务状态查询 URL：<span class="mono">${escapeHtml(statusUrl)}</span>（请保存）</div>
    `;

    const tableBody = root.querySelector('#results-body')!;
    tableBody.innerHTML = renderTableBody(result.items);

    const orderByEl = root.querySelector('#order-by') as HTMLSelectElement;
    if (orderByEl) orderByEl.value = result.order_by;

    const offsetEl = root.querySelector('#offset') as HTMLSpanElement;
    const pageInfoEl = root.querySelector('#page-info') as HTMLSpanElement;
    if (offsetEl) offsetEl.textContent = `${result.offset}`;
    if (pageInfoEl) {
        const start = result.total === 0 ? 0 : result.offset + 1;
        const end = Math.min(result.total, result.offset + result.limit);
        pageInfoEl.textContent = `${start}-${end} / ${result.total}`;
    }

    const prevBtn = root.querySelector('#btn-prev') as HTMLButtonElement;
    const nextBtn = root.querySelector('#btn-next') as HTMLButtonElement;
    prevBtn.disabled = result.offset <= 0;
    nextBtn.disabled = result.offset + result.limit >= result.total;
}

function renderBase(root: HTMLElement, jobId: string, currentUrl: string) {
    root.innerHTML = `
        <div class="status-row">
            <span>Job ID:</span>
            <span class="mono">${escapeHtml(jobId)}</span>
            <span>Status:</span>
            <span class="status-badge queued" id="status-badge">loading</span>
            <button class="btn-secondary" id="btn-refresh" type="button">刷新</button>
            <button class="btn-secondary" id="btn-copy" type="button">复制 URL</button>
        </div>

        <div class="notice">
            结果页面 URL（可保存/分享）：<span class="mono">${escapeHtml(currentUrl)}</span>
        </div>

        <div class="form-grid" style="margin-top:12px;">
            <div class="form-row">
                <label for="order-by">排序（order_by）</label>
                <div style="display:flex; gap:10px; align-items:center; flex-wrap:wrap;">
                    <select id="order-by" class="select" style="max-width: 260px;">
                        <option value="fdr_ranksum">fdr_ranksum (default)</option>
                        <option value="fdr_ttest">fdr_ttest</option>
                        <option value="p_ranksum">p_ranksum</option>
                        <option value="p_ttest">p_ttest</option>
                        <option value="abs_log2fc">abs_log2fc</option>
                    </select>
                    <span class="form-help">offset=<span id="offset" class="mono">0</span></span>
                </div>
            </div>

            <div class="form-row">
                <label>分页</label>
                <div style="display:flex; gap:10px; align-items:center; flex-wrap:wrap;">
                    <input id="limit" class="input" type="number" min="1" max="500" value="50" style="max-width: 120px;" />
                    <button class="btn-secondary" id="btn-apply" type="button">应用</button>
                    <button class="btn-secondary" id="btn-prev" type="button">上一页</button>
                    <button class="btn-secondary" id="btn-next" type="button">下一页</button>
                    <span id="page-info" class="form-help"></span>
                </div>
            </div>
        </div>

        <div id="summary" class="notice"></div>
        <div id="error" style="margin-top:12px;"></div>

        <div class="results-table-container" style="margin-top:12px;">
            <table class="results-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>log2fc</th>
                        <th>mean1</th>
                        <th>mean2</th>
                        <th>pct1</th>
                        <th>pct2</th>
                        <th>p_ttest</th>
                        <th>fdr_ttest</th>
                        <th>p_ranksum</th>
                        <th>fdr_ranksum</th>
                    </tr>
                </thead>
                <tbody id="results-body"></tbody>
            </table>
        </div>
    `;

    const copyBtn = root.querySelector('#btn-copy') as HTMLButtonElement | null;
    copyBtn?.addEventListener('click', async () => {
        try {
            await navigator.clipboard.writeText(window.location.href);
            copyBtn.textContent = '已复制';
            window.setTimeout(() => (copyBtn.textContent = '复制 URL'), 1200);
        } catch {
            copyBtn.textContent = '复制失败';
            window.setTimeout(() => (copyBtn.textContent = '复制 URL'), 1200);
        }
    });
}

async function init() {
    new ThemeManager();

    const root = document.getElementById('de-result-root');
    if (!root) return;

    const url = new URL(window.location.href);
    const dataset = url.searchParams.get('dataset') || '';
    const jobId = url.searchParams.get('job_id') || '';

    if (!dataset || !jobId) {
        root.innerHTML = `<div class="error-box">缺少参数：需要 <span class="mono">dataset</span> 和 <span class="mono">job_id</span>。</div>`;
        return;
    }

    let datasetsInfo: DatasetsResponse | null = null;
    try {
        datasetsInfo = await fetchDatasets();
    } catch {
        // Optional.
    }

    if (datasetsInfo?.title) {
        const titleEl = document.querySelector('.header-title');
        if (titleEl) titleEl.textContent = datasetsInfo.title;
        document.title = `${datasetsInfo.title} - DE Result`;
    }

    const datasetInfoContainer = document.getElementById('dataset-info');
    if (datasetInfoContainer && datasetsInfo) {
        createDatasetSelector(datasetInfoContainer, datasetsInfo.datasets, dataset);
    }

    const titleLink = document.querySelector('.header-title') as HTMLAnchorElement | null;
    if (titleLink) {
        titleLink.href = `/?dataset=${encodeURIComponent(dataset)}`;
    }

    const backBtn = document.getElementById('btn-back') as HTMLAnchorElement | null;
    if (backBtn) {
        backBtn.href = `/?dataset=${encodeURIComponent(dataset)}`;
    }

    const jobUrl = new URL(window.location.href);
    jobUrl.pathname = '/de_job.html';
    jobUrl.searchParams.set('dataset', dataset);
    jobUrl.searchParams.set('job_id', jobId);
    jobUrl.searchParams.delete('offset');
    jobUrl.searchParams.delete('limit');
    jobUrl.searchParams.delete('order_by');

    const jobBtn = document.getElementById('btn-job') as HTMLAnchorElement | null;
    if (jobBtn) {
        jobBtn.href = jobUrl.toString();
    }

    const api = new ApiClient(`/d/${encodeURIComponent(dataset)}/api`);

    renderBase(root, jobId, window.location.href);

    const badge = root.querySelector('#status-badge') as HTMLElement;
    const refreshBtn = root.querySelector('#btn-refresh') as HTMLButtonElement;
    const orderByEl = root.querySelector('#order-by') as HTMLSelectElement;
    const limitEl = root.querySelector('#limit') as HTMLInputElement;
    const applyBtn = root.querySelector('#btn-apply') as HTMLButtonElement;
    const prevBtn = root.querySelector('#btn-prev') as HTMLButtonElement;
    const nextBtn = root.querySelector('#btn-next') as HTMLButtonElement;
    const errorEl = root.querySelector('#error') as HTMLDivElement;

    const getQueryState = () => {
        const u = new URL(window.location.href);
        const offset = Math.max(0, parseInt(u.searchParams.get('offset') || '0', 10) || 0);
        const limit = Math.max(1, Math.min(500, parseInt(u.searchParams.get('limit') || limitEl.value || '50', 10) || 50));
        const order_by = (u.searchParams.get('order_by') as SomaDeResultOrderBy) || (orderByEl.value as SomaDeResultOrderBy) || 'fdr_ranksum';
        return { offset, limit, order_by };
    };

    const setQueryState = (next: { offset?: number; limit?: number; order_by?: SomaDeResultOrderBy }) => {
        const u = new URL(window.location.href);
        if (next.offset !== undefined) u.searchParams.set('offset', String(next.offset));
        if (next.limit !== undefined) u.searchParams.set('limit', String(next.limit));
        if (next.order_by !== undefined) u.searchParams.set('order_by', next.order_by);
        window.history.replaceState({}, '', u.toString());
    };

    let lastStatus: SomaDeJobStatusResponse | null = null;
    let lastResult: SomaDeJobResultResponse | null = null;
    let inFlight = false;

    const renderError = (message: string | null) => {
        if (!message) {
            errorEl.innerHTML = '';
            return;
        }
        errorEl.innerHTML = `<div class="error-box">${escapeHtml(message)}</div>`;
    };

    const refreshStatus = async (): Promise<SomaDeJobStatusResponse> => {
        const info = await api.getSomaDeJob(jobId);
        badge.textContent = info.status;
        setBadgeClass(badge, info.status);
        lastStatus = info;
        return info;
    };

    const loadResults = async () => {
        const { offset, limit, order_by } = getQueryState();
        setQueryState({ offset, limit, order_by });
        limitEl.value = String(limit);
        orderByEl.value = order_by;

        const statusUrl = jobUrl.toString();

        const status = lastStatus || (await refreshStatus());
        if (status.status !== 'completed') {
            renderError(`任务未完成（status: ${status.status}）。请前往“任务进度”页面等待完成。`);
            return;
        }

        const result = await api.getSomaDeJobResult(jobId, { offset, limit, order_by });
        lastResult = result;
        renderError(null);
        renderResult(root, result, statusUrl);
    };

    const refreshAll = async () => {
        if (inFlight) return;
        inFlight = true;
        refreshBtn.disabled = true;
        try {
            await refreshStatus();
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
        const order_by = orderByEl.value as SomaDeResultOrderBy;
        setQueryState({ offset: 0, limit, order_by });
        void refreshAll();
    });

    orderByEl.addEventListener('change', () => {
        const limit = Math.max(1, Math.min(500, parseInt(limitEl.value || '50', 10) || 50));
        setQueryState({ offset: 0, limit, order_by: orderByEl.value as SomaDeResultOrderBy });
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

