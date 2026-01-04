// Cell Query Panel Component
// Allows querying bins/cells expressing a gene

import { ApiClient, BinExpressionInfo, BinQueryResult, GeneCategoryMeanItem, GeneStats } from '../api/client';

export interface CellQueryPanelConfig {
    onBinSelect?: (bin: BinExpressionInfo) => void;
}

export class CellQueryPanel {
    private container: HTMLElement;
    private api: ApiClient;
    private config: CellQueryPanelConfig;

    private currentGene: string | null = null;
    private currentCategoryColumn: string | null = null;
    private currentThreshold: number = 0;
    private currentResult: BinQueryResult | null = null;
    private stats: GeneStats | null = null;
    private categoryMeans: GeneCategoryMeanItem[] | null = null;
    private categoryMeansRequestId: number = 0;

    constructor(
        container: HTMLElement,
        api: ApiClient,
        config: CellQueryPanelConfig = {}
    ) {
        this.container = container;
        this.api = api;
        this.config = config;
        this.render();
    }

    private render(): void {
        this.container.innerHTML = `
            <div class="cell-query-panel">
                <div class="query-controls">
                    <div class="threshold-control">
                        <label>Threshold:</label>
                        <input type="range"
                               id="threshold-slider"
                               min="0" max="10" step="0.1" value="0">
                        <span id="threshold-value">0</span>
                    </div>
                    <button id="btn-query" class="btn-primary" disabled>
                        Query Bins
                    </button>
                </div>

                <div id="query-stats" class="query-stats hidden">
                    <div class="stat-row">
                        <span class="stat-label">Expressing:</span>
                        <span id="stat-expressing">-</span>
                    </div>
                    <div class="stat-row">
                        <span class="stat-label">Mean:</span>
                        <span id="stat-mean">-</span>
                    </div>
                    <div class="stat-row">
                        <span class="stat-label">Max:</span>
                        <span id="stat-max">-</span>
                    </div>
                </div>

                <div id="category-mean-hist" class="category-mean-hist hidden">
                    <div class="category-mean-header">
                        <span class="category-mean-title">By category</span>
                        <span id="category-mean-column" class="category-mean-column"></span>
                    </div>
                    <div class="category-mean-chart">
                        <div id="category-mean-bars" class="category-mean-bars"></div>
                    </div>
                    <div id="category-mean-footer" class="category-mean-footer"></div>
                </div>

                <div id="query-results" class="query-results hidden">
                    <div class="results-header">
                        <span id="results-count">0 bins</span>
                    </div>
                    <div id="bin-table-container" class="bin-table-container">
                        <table class="bin-table">
                            <thead>
                                <tr>
                                    <th>Expression</th>
                                    <th>Cells</th>
                                    <th>Position</th>
                                </tr>
                            </thead>
                            <tbody id="bin-table-body"></tbody>
                        </table>
                    </div>
                    <div class="pagination">
                        <button id="btn-prev" disabled>&lt;</button>
                        <span id="page-info">1</span>
                        <button id="btn-next" disabled>&gt;</button>
                    </div>
                </div>
            </div>
        `;

        this.setupEvents();
    }

    private setupEvents(): void {
        // Threshold slider
        const slider = this.container.querySelector('#threshold-slider') as HTMLInputElement;
        const valueDisplay = this.container.querySelector('#threshold-value')!;
        slider.addEventListener('input', () => {
            this.currentThreshold = parseFloat(slider.value);
            valueDisplay.textContent = this.currentThreshold.toFixed(1);
        });

        // Query button
        const queryBtn = this.container.querySelector('#btn-query') as HTMLButtonElement;
        queryBtn.addEventListener('click', () => this.executeQuery());

        // Pagination
        const prevBtn = this.container.querySelector('#btn-prev') as HTMLButtonElement;
        const nextBtn = this.container.querySelector('#btn-next') as HTMLButtonElement;
        prevBtn.addEventListener('click', () => this.loadPage(-1));
        nextBtn.addEventListener('click', () => this.loadPage(1));
    }

    /**
     * Set the gene to query
     */
    async setGene(gene: string): Promise<void> {
        this.currentGene = gene;
        this.currentResult = null;
        this.categoryMeans = null;

        // Enable query button
        const queryBtn = this.container.querySelector('#btn-query') as HTMLButtonElement;
        queryBtn.disabled = false;

        // Load gene statistics
        try {
            this.stats = await this.api.getGeneStats(gene);
            this.showStats();

            // Adjust slider range based on max expression
            const slider = this.container.querySelector('#threshold-slider') as HTMLInputElement;
            slider.max = Math.ceil(this.stats.max_expression).toString();
        } catch (error) {
            console.error('Failed to load gene stats:', error);
        }

        await this.loadCategoryMeans();
    }

    /**
     * Set the category column used for the per-category mean plot.
     */
    setCategoryColumn(column: string): void {
        this.currentCategoryColumn = column;
        void this.loadCategoryMeans();
    }

    private showStats(): void {
        if (!this.stats) return;

        const statsDiv = this.container.querySelector('#query-stats')!;
        statsDiv.classList.remove('hidden');

        this.container.querySelector('#stat-expressing')!.textContent =
            `${this.stats.expressing_bins.toLocaleString()} / ${this.stats.total_bins.toLocaleString()} bins`;
        this.container.querySelector('#stat-mean')!.textContent =
            this.stats.mean_expression.toFixed(3);
        this.container.querySelector('#stat-max')!.textContent =
            this.stats.max_expression.toFixed(3);
    }

    private async loadCategoryMeans(): Promise<void> {
        const histDiv = this.container.querySelector('#category-mean-hist') as HTMLElement | null;
        if (!histDiv) return;

        if (!this.currentGene || !this.currentCategoryColumn) {
            histDiv.classList.add('hidden');
            return;
        }

        const requestId = ++this.categoryMeansRequestId;

        histDiv.classList.remove('hidden');
        const columnEl = this.container.querySelector('#category-mean-column') as HTMLElement | null;
        if (columnEl) columnEl.textContent = this.currentCategoryColumn;

        const barsEl = this.container.querySelector('#category-mean-bars') as HTMLElement | null;
        const footerEl = this.container.querySelector('#category-mean-footer') as HTMLElement | null;
        if (barsEl) barsEl.innerHTML = `<div class="category-mean-loading">Loading…</div>`;
        if (footerEl) footerEl.textContent = '';

        try {
            const items = await this.api.getGeneCategoryMeans(this.currentGene, this.currentCategoryColumn);
            if (requestId !== this.categoryMeansRequestId) return;

            this.categoryMeans = items;
            this.showCategoryMeans();
        } catch (error) {
            if (requestId !== this.categoryMeansRequestId) return;
            console.error('Failed to load category means:', error);
            if (barsEl) barsEl.innerHTML = `<div class="category-mean-error">Failed to load category means</div>`;
        }
    }

    private showCategoryMeans(): void {
        const histDiv = this.container.querySelector('#category-mean-hist') as HTMLElement | null;
        const barsEl = this.container.querySelector('#category-mean-bars') as HTMLElement | null;
        const footerEl = this.container.querySelector('#category-mean-footer') as HTMLElement | null;
        if (!histDiv || !barsEl || !footerEl || !this.categoryMeans) return;

        if (this.categoryMeans.length === 0) {
            barsEl.innerHTML = `<div class="category-mean-empty">No categories</div>`;
            footerEl.textContent = '';
            return;
        }

        const maxMean = Math.max(0, ...this.categoryMeans.map(d => d.mean_expression));
        const scaleMax = maxMean > 0 ? maxMean : 1;

        barsEl.innerHTML = '';
        for (const d of this.categoryMeans) {
            const heightPct = (d.mean_expression / scaleMax) * 100;

            const row = document.createElement('div');
            row.className = 'category-mean-row';
            row.title = `${d.value}\nmean: ${d.mean_expression.toFixed(4)}\nbins: ${d.bin_count.toLocaleString()}`;

            const track = document.createElement('div');
            track.className = 'category-mean-row-bar-track';

            const fill = document.createElement('div');
            fill.className = 'category-mean-row-bar-fill';
            fill.style.height = `${heightPct.toFixed(2)}%`;
            fill.style.backgroundColor = d.color || 'var(--accent-secondary)';

            const label = document.createElement('div');
            label.className = 'category-mean-row-label';
            label.textContent = d.value;

            const value = document.createElement('div');
            value.className = 'category-mean-row-value';
            value.textContent = d.mean_expression.toFixed(4);

            track.appendChild(fill);
            row.append(track, label, value);
            barsEl.appendChild(row);
        }

        footerEl.textContent = `max mean: ${maxMean.toFixed(4)}`;
        histDiv.classList.remove('hidden');
    }

    private async executeQuery(offset: number = 0): Promise<void> {
        if (!this.currentGene) return;

        try {
            this.currentResult = await this.api.getBinsExpressingGene(
                this.currentGene,
                {
                    threshold: this.currentThreshold,
                    offset,
                    limit: 50,
                }
            );

            this.showResults();
        } catch (error) {
            console.error('Query failed:', error);
        }
    }

    private showResults(): void {
        if (!this.currentResult) return;

        const resultsDiv = this.container.querySelector('#query-results')!;
        resultsDiv.classList.remove('hidden');

        // Update count
        this.container.querySelector('#results-count')!.textContent =
            `${this.currentResult.total_count.toLocaleString()} bins`;

        // Update table
        const tbody = this.container.querySelector('#bin-table-body')!;
        tbody.innerHTML = this.currentResult.bins.map(bin => `
            <tr data-bin-index="${bin.bin_index}" data-x="${bin.bin_x}" data-y="${bin.bin_y}">
                <td>${bin.expression.toFixed(3)}</td>
                <td>${bin.cell_count}</td>
                <td>(${bin.bin_x}, ${bin.bin_y})</td>
            </tr>
        `).join('');

        // Add row click events
        tbody.querySelectorAll('tr').forEach(row => {
            row.addEventListener('click', () => {
                const binIndex = parseInt((row as HTMLElement).dataset.binIndex!);
                const bin = this.currentResult!.bins.find(b => b.bin_index === binIndex);
                if (bin && this.config.onBinSelect) {
                    this.config.onBinSelect(bin);
                }
            });
        });

        // Update pagination
        this.updatePagination();
    }

    private updatePagination(): void {
        if (!this.currentResult) return;

        const { offset, limit, total_count } = this.currentResult;
        const currentPage = Math.floor(offset / limit) + 1;
        const totalPages = Math.ceil(total_count / limit);

        const prevBtn = this.container.querySelector('#btn-prev') as HTMLButtonElement;
        const nextBtn = this.container.querySelector('#btn-next') as HTMLButtonElement;
        const pageInfo = this.container.querySelector('#page-info')!;

        prevBtn.disabled = currentPage <= 1;
        nextBtn.disabled = currentPage >= totalPages;
        pageInfo.textContent = `${currentPage}/${totalPages}`;
    }

    private loadPage(direction: number): void {
        if (!this.currentResult) return;

        const newOffset = this.currentResult.offset + (direction * this.currentResult.limit);
        if (newOffset >= 0 && newOffset < this.currentResult.total_count) {
            this.executeQuery(newOffset);
        }
    }

    /**
     * Clear the query panel
     */
    clear(): void {
        this.currentGene = null;
        this.currentResult = null;
        this.stats = null;
        this.categoryMeans = null;
        this.categoryMeansRequestId++;

        const queryBtn = this.container.querySelector('#btn-query') as HTMLButtonElement;
        queryBtn.disabled = true;

        this.container.querySelector('#query-stats')!.classList.add('hidden');
        this.container.querySelector('#category-mean-hist')!.classList.add('hidden');
        this.container.querySelector('#query-results')!.classList.add('hidden');
    }
}
