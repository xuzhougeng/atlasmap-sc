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
import { VolcanoPlot } from './components/VolcanoPlot';

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
        return `<tr><td colspan="10" style="color:var(--text-secondary); padding:16px; text-align:center;">No results</td></tr>`;
    }

    return items
        .map((it) => {
            return `
                <tr>
                    <td>${escapeHtml(it.gene)}</td>
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

type VolcanoMetric = 'fdr_ranksum' | 'p_ranksum' | 'fdr_ttest' | 'p_ttest';

function getMetricValue(item: SomaDeResultItem, metric: VolcanoMetric): number {
    switch (metric) {
        case 'fdr_ranksum':
            return item.fdr_ranksum;
        case 'p_ranksum':
            return item.p_ranksum;
        case 'fdr_ttest':
            return item.fdr_ttest;
        case 'p_ttest':
            return item.p_ttest;
    }
}

function negLog10(v: number): number {
    const eps = 1e-300;
    const safe = Math.max(eps, Math.min(1, v));
    return -Math.log10(safe);
}

function renderResult(root: HTMLElement, result: SomaDeJobResultResponse) {
    const params = result.params;
    const group2Label = params.group2?.length ? params.group2.join(', ') : '(one-vs-rest)';

    // Compact Summary Strip
    root.querySelector('#summary')!.innerHTML = `
        <div class="de-result-summary-item">
            <span class="label">Group By</span>
            <span class="value">${escapeHtml(params.groupby)}</span>
        </div>
        <div class="de-result-summary-item">
            <span class="label">Group 1</span>
            <span class="value" title="${escapeHtml(params.group1.join(', '))}">${escapeHtml(params.group1.join(', '))}</span>
        </div>
        <div class="de-result-summary-item">
            <span class="label">Group 2</span>
            <span class="value" title="${escapeHtml(group2Label)}">${escapeHtml(group2Label)}</span>
        </div>
        <div class="de-result-summary-divider"></div>
        <div class="de-result-summary-item">
            <span class="label">Tests</span>
            <span class="value">${escapeHtml((params.tests || []).join(', '))}</span>
        </div>
        <div class="de-result-summary-item">
            <span class="label">N1 / N2</span>
            <span class="value">${result.n1} / ${result.n2}</span>
        </div>
        <div class="de-result-summary-item">
            <span class="label">Total Genes</span>
            <span class="value">${result.total.toLocaleString()}</span>
        </div>
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

function renderBase(root: HTMLElement, jobId: string) {
    // Add result page class to body
    document.body.classList.add('de-result-page');

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

        <!-- Main Split Panel -->
        <div class="de-result-main">
            <!-- Volcano Plot Panel -->
            <div class="de-result-panel de-result-volcano-panel">
                <div class="de-result-panel-header">
                    <span class="de-result-panel-title">Volcano Plot</span>
                    <div class="de-result-panel-controls">
                        <div class="de-result-control-group">
                            <span class="label">Metric</span>
                            <select id="plot-metric" class="de-result-select">
                                <option value="fdr_ranksum">fdr_ranksum</option>
                                <option value="p_ranksum">p_ranksum</option>
                                <option value="fdr_ttest">fdr_ttest</option>
                                <option value="p_ttest">p_ttest</option>
                            </select>
                        </div>
                        <div class="de-result-control-group">
                            <span class="label">|log2FC| ≥</span>
                            <input id="plot-fc" class="de-result-input" type="number" step="0.1" value="1" />
                        </div>
	                        <div class="de-result-control-group">
	                            <span class="label">Sig ≤</span>
	                            <input id="plot-sig" class="de-result-input" type="number" step="0.001" min="0" max="1" value="0.05" />
	                        </div>
	                        <div class="de-result-control-group">
	                            <input id="plot-exclude-zero" type="checkbox" checked />
	                            <span class="label">Hide pct1=pct2=0</span>
	                        </div>
	                        <button class="de-result-btn-sm primary" id="btn-load-plot" type="button">Load Plot</button>
	                        <span id="plot-status" class="de-result-panel-status"></span>
	                    </div>
	                </div>
                <div class="de-result-volcano-body">
                    <canvas id="volcano-gl"></canvas>
                    <canvas id="volcano-overlay"></canvas>
                    <div class="de-result-volcano-placeholder" id="volcano-placeholder">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                            <path d="M12 2L2 22h20L12 2z"/>
                            <circle cx="12" cy="16" r="1"/>
                            <line x1="12" y1="9" x2="12" y2="13"/>
                        </svg>
                        <p>Click "Load Plot" to fetch all genes and render the volcano plot.</p>
                        <p class="hint">(This may take a moment for large datasets)</p>
                    </div>
                </div>
                <div id="plot-error" class="de-result-volcano-error" style="display:none;"></div>
            </div>

            <!-- Table Panel -->
            <div class="de-result-panel de-result-table-panel">
                <div class="de-result-panel-header">
                    <span class="de-result-panel-title">Gene Results</span>
                    <div class="de-result-panel-controls">
                        <div class="de-result-control-group">
                            <span class="label">Sort By</span>
                            <select id="order-by" class="de-result-select">
                                <option value="fdr_ranksum">fdr_ranksum</option>
                                <option value="fdr_ttest">fdr_ttest</option>
                                <option value="p_ranksum">p_ranksum</option>
                                <option value="p_ttest">p_ttest</option>
                                <option value="abs_log2fc">abs_log2fc</option>
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
                    <table class="de-result-table">
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
                <div class="de-result-table-footer">
                    <div class="de-result-pagination">
                        <button class="de-result-btn-sm" id="btn-prev" type="button">Prev</button>
                        <span id="page-info" class="de-result-pagination-info"></span>
                        <button class="de-result-btn-sm" id="btn-next" type="button">Next</button>
                    </div>
                </div>
            </div>
        </div>
    `;

    const copyBtn = root.querySelector('#btn-copy') as HTMLButtonElement | null;
    copyBtn?.addEventListener('click', async () => {
        try {
            await navigator.clipboard.writeText(window.location.href);
            copyBtn.textContent = 'Copied';
            window.setTimeout(() => (copyBtn.textContent = 'Copy Link'), 1200);
        } catch {
            copyBtn.textContent = 'Copy failed';
            window.setTimeout(() => (copyBtn.textContent = 'Copy Link'), 1200);
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
        root.innerHTML = `<div class="error-box">Missing parameters: requires <span class="mono">dataset</span> and <span class="mono">job_id</span>.</div>`;
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

    renderBase(root, jobId);

    const badge = root.querySelector('#status-badge') as HTMLElement;
    const refreshBtn = root.querySelector('#btn-refresh') as HTMLButtonElement;
    const orderByEl = root.querySelector('#order-by') as HTMLSelectElement;
    const limitEl = root.querySelector('#limit') as HTMLInputElement;
    const applyBtn = root.querySelector('#btn-apply') as HTMLButtonElement;
    const prevBtn = root.querySelector('#btn-prev') as HTMLButtonElement;
    const nextBtn = root.querySelector('#btn-next') as HTMLButtonElement;
    const errorEl = root.querySelector('#error') as HTMLDivElement;

    const plotMetricEl = root.querySelector('#plot-metric') as HTMLSelectElement;
    const plotFcEl = root.querySelector('#plot-fc') as HTMLInputElement;
    const plotSigEl = root.querySelector('#plot-sig') as HTMLInputElement;
    const plotExcludeZeroEl = root.querySelector('#plot-exclude-zero') as HTMLInputElement;
    const plotBtn = root.querySelector('#btn-load-plot') as HTMLButtonElement;
    const plotStatusEl = root.querySelector('#plot-status') as HTMLSpanElement;
    const plotPlaceholderEl = root.querySelector('#volcano-placeholder') as HTMLDivElement;
    const plotErrorEl = root.querySelector('#plot-error') as HTMLDivElement;
    const plotGlCanvas = root.querySelector('#volcano-gl') as HTMLCanvasElement;
    const plotOverlayCanvas = root.querySelector('#volcano-overlay') as HTMLCanvasElement;

    let volcanoPlot: VolcanoPlot | null = null;
    let plotItems: SomaDeResultItem[] | null = null;
    let plotLoading = false;
    let plotSupported = true;

    const renderPlotError = (message: string | null) => {
        if (!message) {
            plotErrorEl.innerHTML = '';
            plotErrorEl.style.display = 'none';
            return;
        }
        plotErrorEl.style.display = 'block';
        plotErrorEl.innerHTML = `<div class="de-result-error-box">${escapeHtml(message)}</div>`;
    };

    const getPlotThresholds = () => {
        const fc = Math.max(0, parseFloat(plotFcEl.value || '1') || 1);
        const sigCutoffRaw = parseFloat(plotSigEl.value || '0.05');
        const sigCutoff = Number.isFinite(sigCutoffRaw) ? sigCutoffRaw : 0.05;
        const sigY = negLog10(sigCutoff);
        return { fc, sigCutoff, sigY };
    };

    const ensurePlotInstance = () => {
        if (!plotSupported) return false;
        if (volcanoPlot) return true;
        try {
            volcanoPlot = new VolcanoPlot(plotGlCanvas, plotOverlayCanvas, { pointSize: 3 });
            const metric = plotMetricEl.value as VolcanoMetric;
            volcanoPlot.setLabels('log2FC', `-log10(${metric})`);
            const { fc, sigY } = getPlotThresholds();
            volcanoPlot.setThresholds(fc, sigY);
            volcanoPlot.setData([]);
            return true;
        } catch (e: any) {
            console.error(e);
            renderPlotError(e?.message || String(e));
            plotSupported = false;
            plotBtn.disabled = true;
            plotPlaceholderEl.textContent = 'WebGL is not available; volcano plot cannot be rendered.';
            return false;
        }
    };

    const updatePlotThresholds = () => {
        if (!volcanoPlot) return;
        const { fc, sigY } = getPlotThresholds();
        volcanoPlot.setThresholds(fc, sigY);
    };

    const updatePlotData = () => {
        if (!plotItems) return;
        if (!ensurePlotInstance() || !volcanoPlot) return;

        const excludeZero = plotExcludeZeroEl.checked;
        const filtered = excludeZero
            ? plotItems.filter(it => it.pct1 > 0 || it.pct2 > 0)
            : plotItems;

        const metric = plotMetricEl.value as VolcanoMetric;
        volcanoPlot.setLabels('log2FC', `-log10(${metric})`);

        const points = filtered.map(it => ({
            x: it.log2fc,
            y: negLog10(getMetricValue(it, metric)),
        }));
        volcanoPlot.setData(points);
        updatePlotThresholds();

        plotPlaceholderEl.style.display = 'none';
        plotStatusEl.textContent = `Loaded ${plotItems.length.toLocaleString()} genes (plotted ${filtered.length.toLocaleString()})`;
    };

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
            errorEl.style.display = 'none';
            return;
        }
        errorEl.style.display = 'block';
        errorEl.innerHTML = `<div class="de-result-error-box">${escapeHtml(message)}</div>`;
    };

    const refreshStatus = async (): Promise<SomaDeJobStatusResponse> => {
        const info = await api.getSomaDeJob(jobId);
        badge.textContent = info.status;
        setBadgeClass(badge, info.status);
        lastStatus = info;

        plotBtn.disabled = !plotSupported || plotLoading || info.status !== 'completed';
        const plottedCount = plotItems
            ? (plotExcludeZeroEl.checked
                ? plotItems.filter(it => it.pct1 > 0 || it.pct2 > 0).length
                : plotItems.length)
            : 0;
        plotStatusEl.textContent = !plotSupported
            ? 'WebGL unavailable'
            : info.status === 'completed'
                ? (plotItems ? `Loaded ${plotItems.length.toLocaleString()} genes (plotted ${plottedCount.toLocaleString()})` : 'Ready to load')
                : `Waiting (status: ${info.status})`;

        return info;
    };

    const loadAllResultsForPlot = async (): Promise<void> => {
        if (plotLoading) return;
        const status = lastStatus || (await refreshStatus());
        if (status.status !== 'completed') {
            renderPlotError(`Job is not completed yet (status: ${status.status}).`);
            return;
        }

        if (!ensurePlotInstance()) return;

        plotLoading = true;
        plotBtn.disabled = true;
        renderPlotError(null);
        plotPlaceholderEl.style.display = 'flex';
        plotPlaceholderEl.textContent = 'Loading plot data...';
        plotStatusEl.textContent = 'Fetching...';
        plotBtn.textContent = 'Loading...';

        try {
            const items: SomaDeResultItem[] = [];
            let total = 0;
            let offset = 0;
            const limit = 500;

            while (true) {
                const res = await api.getSomaDeJobResult(jobId, { offset, limit, order_by: 'fdr_ranksum' });
                if (total === 0) total = res.total;
                items.push(...res.items);
                offset += res.items.length;
                if (total > 0) {
                    plotStatusEl.textContent = `Loaded ${Math.min(offset, total).toLocaleString()}/${total.toLocaleString()}`;
                } else {
                    plotStatusEl.textContent = `Loaded ${offset.toLocaleString()}`;
                }
                if (res.items.length === 0 || (total > 0 && offset >= total)) break;
            }

            plotItems = items;
            plotStatusEl.textContent = `Loaded ${items.length.toLocaleString()} genes`;
            plotBtn.textContent = 'Reload plot';
            plotBtn.disabled = false;

            updatePlotData();
        } catch (e: any) {
            console.error(e);
            renderPlotError(e?.message || String(e));
            plotBtn.disabled = false;
            plotBtn.textContent = 'Load plot (all genes)';
            plotStatusEl.textContent = 'Failed to load';
        } finally {
            plotLoading = false;
        }
    };

    const loadResults = async () => {
        const { offset, limit, order_by } = getQueryState();
        setQueryState({ offset, limit, order_by });
        limitEl.value = String(limit);
        orderByEl.value = order_by;

        const status = lastStatus || (await refreshStatus());
        if (status.status !== 'completed') {
            renderError(`Job is not completed yet (status: ${status.status}). Please use the Job Status page and wait for completion.`);
            return;
        }

        const result = await api.getSomaDeJobResult(jobId, { offset, limit, order_by });
        lastResult = result;
        renderError(null);
        renderResult(root, result);
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

    plotBtn.addEventListener('click', () => void loadAllResultsForPlot());

    plotMetricEl.addEventListener('change', () => {
        renderPlotError(null);
        updatePlotData();
    });
    plotFcEl.addEventListener('change', () => {
        renderPlotError(null);
        updatePlotThresholds();
    });
    plotSigEl.addEventListener('change', () => {
        renderPlotError(null);
        updatePlotThresholds();
    });
    plotExcludeZeroEl.addEventListener('change', () => {
        renderPlotError(null);
        updatePlotData();
    });

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
