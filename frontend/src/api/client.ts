// API Client for AtlasMap backend

export interface Metadata {
    dataset_name: string;
    n_cells: number;
    n_genes_total: number;
    n_genes_preaggregated: number;
    zoom_levels: number;
    tile_size: number;
    preaggregated_genes: string[];
    categories: Record<string, CategoryInfo>;
    bounds: {
        min_x: number;
        max_x: number;
        min_y: number;
        max_y: number;
    };
    umap_key?: string;
    coordinate_systems?: CoordinateSystemInfo[];
    default_coordinate_system?: string;
    coordinate_system?: string;
}

export interface CoordinateSystemInfo {
    key: string;
    zarr_path: string;
}

export interface CategoryInfo {
    values: string[];
    mapping: Record<string, number>;
}

export interface GeneInfo {
    name: string;
    index: number;
    preaggregated: boolean;
}

export interface CategoryLegendItem {
    value: string;
    color: string;
    index: number;
    cell_count: number;
}

export interface CategoryCentroidItem {
    value: string;
    color: string;
    index: number;
    bin_count: number;
    cell_count: number;
    x: number | null;
    y: number | null;
}

export interface BinExpressionInfo {
    bin_index: number;
    bin_x: number;
    bin_y: number;
    cell_count: number;
    expression: number;
    categories: Record<string, string>;
}

export interface BinQueryResult {
    bins: BinExpressionInfo[];
    total_count: number;
    offset: number;
    limit: number;
    gene: string;
    threshold: number;
}

export interface GeneStats {
    gene: string;
    index: number;
    expressing_bins: number;
    total_bins: number;
    total_cells: number;
    mean_expression: number;
    max_expression: number;
    p80_expression?: number;
}

export interface GeneCategoryMeanItem {
    value: string;
    color: string;
    index: number;
    bin_count: number;
    mean_expression: number;
}

export interface GeneCategoryMeansResponse {
    gene: string;
    column: string;
    items: GeneCategoryMeanItem[];
}

export interface CellInfo {
    joinid: number;
    x: number;
    y: number;
    expression?: number;
    category?: string;
}

export interface CellQueryResult {
    cells: CellInfo[];
    total_count: number;
    truncated: boolean;
}

export interface BinQueryParams {
    threshold?: number;
    offset?: number;
    limit?: number;
}

export type SomaDeTest = 'ttest' | 'ranksum';
export type SomaDeJobStatus = 'queued' | 'running' | 'completed' | 'failed' | 'cancelled';
export type SomaDeResultOrderBy = 'fdr_ranksum' | 'fdr_ttest' | 'p_ranksum' | 'p_ttest' | 'abs_log2fc';

export interface SomaDeJobCreateRequest {
    groupby: string;
    group1: string[];
    group2: string[];
    tests?: SomaDeTest[];
    max_cells_per_group?: number;
    seed?: number;
    limit?: number;
}

export interface SomaDeJobCreateResponse {
    job_id: string;
    status: SomaDeJobStatus;
}

export interface SomaDeJobStatusResponse {
    job_id: string;
    status: SomaDeJobStatus;
    created_at?: string;
    started_at?: string;
    finished_at?: string;
    progress?: {
        phase?: string;
        done?: number;
        total?: number;
    };
    n1?: number;
    n2?: number;
    error?: string;
}

interface SomaObsColumnsResponse {
    columns: string[];
}

interface SomaObsValuesResponse {
    column: string;
    values: string[];
}

export interface SomaDeJobParams {
    dataset_id: string;
    groupby: string;
    group1: string[];
    group2: string[];
    tests: string[];
    max_cells_per_group: number;
    seed: number;
    limit: number;
}

export interface SomaDeResultItem {
    gene: string;
    gene_joinid: number;
    mean1: number;
    mean2: number;
    pct1: number;
    pct2: number;
    log2fc: number;
    p_ttest: number;
    fdr_ttest: number;
    p_ranksum: number;
    fdr_ranksum: number;
}

export interface SomaDeJobResultResponse {
    params: SomaDeJobParams;
    n1: number;
    n2: number;
    total: number;
    offset: number;
    limit: number;
    order_by: SomaDeResultOrderBy;
    items: SomaDeResultItem[];
}

export class ApiClient {
    private baseUrl: string;
    private defaultQuery: Record<string, string>;

    constructor(baseUrl: string, defaultQuery: Record<string, string> = {}) {
        this.baseUrl = baseUrl;
        this.defaultQuery = defaultQuery;
    }

    private buildUrl(path: string, query?: URLSearchParams): string {
        const base = this.baseUrl.replace(/\/$/, '');
        const merged = new URLSearchParams();

        for (const [k, v] of Object.entries(this.defaultQuery)) {
            if (typeof v !== 'undefined' && v !== null && String(v) !== '') {
                merged.set(k, String(v));
            }
        }

        if (query) {
            for (const [k, v] of query.entries()) {
                merged.set(k, v);
            }
        }

        const suffix = merged.toString() ? `?${merged.toString()}` : '';
        return `${base}${path}${suffix}`;
    }

    private async readError(response: Response): Promise<string> {
        try {
            const text = await response.text();
            return text || response.statusText;
        } catch {
            return response.statusText;
        }
    }

    private async fetchJson<T>(url: string, init?: RequestInit): Promise<T> {
        const response = await fetch(url, init);
        if (!response.ok) {
            const message = await this.readError(response);
            throw new Error(`HTTP ${response.status}: ${message}`);
        }
        return response.json();
    }

    async getMetadata(): Promise<Metadata> {
        return this.fetchJson(this.buildUrl('/metadata'));
    }

    async getGenes(): Promise<{ genes: string[]; total: number }> {
        return this.fetchJson(this.buildUrl('/genes'));
    }

    async getGeneInfo(gene: string): Promise<GeneInfo> {
        return this.fetchJson(this.buildUrl(`/genes/${encodeURIComponent(gene)}`));
    }

    async getCategories(): Promise<Record<string, CategoryInfo>> {
        return this.fetchJson(this.buildUrl('/categories'));
    }

    async getStats(): Promise<Record<string, any>> {
        return this.fetchJson(this.buildUrl('/stats'));
    }

    /**
     * Search genes by prefix
     */
    async searchGenes(query: string, limit: number = 20): Promise<string[]> {
        // Client-side filtering (since genes list is typically small)
        const { genes } = await this.getGenes();
        const queryLower = query.toLowerCase();
        return genes
            .filter(g => g.toLowerCase().includes(queryLower))
            .slice(0, limit);
    }

    /**
     * Get color mapping for a category column
     */
    async getCategoryColors(column: string): Promise<Record<string, string>> {
        return this.fetchJson(this.buildUrl(`/categories/${encodeURIComponent(column)}/colors`));
    }

    /**
     * Get legend data for a category column
     */
    async getCategoryLegend(column: string): Promise<CategoryLegendItem[]> {
        return this.fetchJson(this.buildUrl(`/categories/${encodeURIComponent(column)}/legend`));
    }

    async getCategoryCentroids(column: string): Promise<CategoryCentroidItem[]> {
        return this.fetchJson(this.buildUrl(`/categories/${encodeURIComponent(column)}/centroids`));
    }

    /**
     * Get available category columns
     */
    getAvailableCategories(metadata: Metadata): string[] {
        return Object.keys(metadata.categories);
    }

    /**
     * Get bins expressing a gene above threshold
     */
    async getBinsExpressingGene(
        gene: string,
        params: BinQueryParams = {}
    ): Promise<BinQueryResult> {
        const query = new URLSearchParams();
        if (params.threshold !== undefined) {
            query.set('threshold', params.threshold.toString());
        }
        if (params.offset !== undefined) {
            query.set('offset', params.offset.toString());
        }
        if (params.limit !== undefined) {
            query.set('limit', params.limit.toString());
        }

        return this.fetchJson(this.buildUrl(`/genes/${encodeURIComponent(gene)}/bins`, query));
    }

    /**
     * Get gene expression statistics
     * @param gene Gene name
     * @param zoom Optional zoom level (defaults to 0 if omitted)
     */
    async getGeneStats(gene: string, zoom?: number): Promise<GeneStats> {
        const query = new URLSearchParams();
        if (zoom !== undefined) {
            query.set('zoom', zoom.toString());
        }
        return this.fetchJson(this.buildUrl(`/genes/${encodeURIComponent(gene)}/stats`, query));
    }

    /**
     * Get mean expression per category value for a gene
     */
    async getGeneCategoryMeans(gene: string, column: string): Promise<GeneCategoryMeanItem[]> {
        const data: GeneCategoryMeansResponse = await this.fetchJson(
            this.buildUrl(`/genes/${encodeURIComponent(gene)}/category/${encodeURIComponent(column)}/means`)
        );
        return data.items;
    }

    async getSomaObsColumns(): Promise<string[]> {
        const data: SomaObsColumnsResponse = await this.fetchJson(this.buildUrl('/soma/obs/columns'));
        return data.columns;
    }

    async getSomaObsColumnValues(column: string): Promise<string[]> {
        const data: SomaObsValuesResponse = await this.fetchJson(
            this.buildUrl(`/soma/obs/${encodeURIComponent(column)}/values`)
        );
        return data.values;
    }

    async createSomaDeJob(params: SomaDeJobCreateRequest): Promise<SomaDeJobCreateResponse> {
        return this.fetchJson(this.buildUrl('/soma/de/jobs'), {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(params),
        });
    }

    async getSomaDeJob(jobId: string): Promise<SomaDeJobStatusResponse> {
        return this.fetchJson(this.buildUrl(`/soma/de/jobs/${encodeURIComponent(jobId)}`));
    }

    async getSomaDeJobResult(
        jobId: string,
        params: { offset?: number; limit?: number; order_by?: SomaDeResultOrderBy } = {}
    ): Promise<SomaDeJobResultResponse> {
        const query = new URLSearchParams();
        if (params.offset !== undefined) query.set('offset', params.offset.toString());
        if (params.limit !== undefined) query.set('limit', params.limit.toString());
        if (params.order_by !== undefined) query.set('order_by', params.order_by);
        return this.fetchJson(this.buildUrl(`/soma/de/jobs/${encodeURIComponent(jobId)}/result`, query));
    }

    async cancelSomaDeJob(jobId: string): Promise<void> {
        const response = await fetch(this.buildUrl(`/soma/de/jobs/${encodeURIComponent(jobId)}`), {
            method: 'DELETE',
        });
        if (!response.ok) {
            const message = await this.readError(response);
            throw new Error(`HTTP ${response.status}: ${message}`);
        }
    }

    /**
     * Get cells within a bounding box for cell-level rendering.
     * @param bounds - The bounding box for the query
     * @param options - Query options including:
     *   - gene: Gene name for expression values
     *   - category: Category column name
     *   - categories: Category filter values (null = no filter, [] = show none)
     *   - limit: Maximum cells to return (default 5000)
     *   - seed: Random seed for deterministic downsampling (default 0)
     * @param signal - AbortSignal for cancellation
     */
    async getSomaCellsInBounds(
        bounds: { minX: number; minY: number; maxX: number; maxY: number },
        options: {
            gene?: string;
            category?: string;
            categories?: string[] | null;
            limit?: number;
            seed?: number;
        } = {},
        signal?: AbortSignal
    ): Promise<CellQueryResult> {
        const query = new URLSearchParams();
        query.set('min_x', bounds.minX.toString());
        query.set('min_y', bounds.minY.toString());
        query.set('max_x', bounds.maxX.toString());
        query.set('max_y', bounds.maxY.toString());

        if (options.gene) {
            query.set('gene', options.gene);
        }
        if (options.category) {
            query.set('category', options.category);
        }
        // Match tile filtering semantics:
        // - undefined => no filter (don't send param)
        // - null      => no filter (explicit, don't send param)
        // - []        => filter-to-none (send "[]")
        // - ["A",..]  => filter to these values
        if (options.categories !== undefined && options.categories !== null) {
            query.set('categories', JSON.stringify(options.categories));
        }
        if (options.limit) {
            query.set('limit', options.limit.toString());
        }
        // Seed for deterministic downsampling (0 = default deterministic sampling)
        if (options.seed !== undefined) {
            query.set('seed', options.seed.toString());
        }

        const url = this.buildUrl('/soma/cells', query);
        const response = await fetch(url, { signal });
        if (!response.ok) {
            const message = await this.readError(response);
            throw new Error(`HTTP ${response.status}: ${message}`);
        }
        return response.json();
    }
}
