// API Client for SOMA-Tiles backend

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

    constructor(baseUrl: string) {
        this.baseUrl = baseUrl;
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
        return this.fetchJson(`${this.baseUrl}/metadata`);
    }

    async getGenes(): Promise<{ genes: string[]; total: number }> {
        return this.fetchJson(`${this.baseUrl}/genes`);
    }

    async getGeneInfo(gene: string): Promise<GeneInfo> {
        return this.fetchJson(`${this.baseUrl}/genes/${encodeURIComponent(gene)}`);
    }

    async getCategories(): Promise<Record<string, CategoryInfo>> {
        return this.fetchJson(`${this.baseUrl}/categories`);
    }

    async getStats(): Promise<Record<string, any>> {
        return this.fetchJson(`${this.baseUrl}/stats`);
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
        return this.fetchJson(`${this.baseUrl}/categories/${encodeURIComponent(column)}/colors`);
    }

    /**
     * Get legend data for a category column
     */
    async getCategoryLegend(column: string): Promise<CategoryLegendItem[]> {
        return this.fetchJson(`${this.baseUrl}/categories/${encodeURIComponent(column)}/legend`);
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

        return this.fetchJson(`${this.baseUrl}/genes/${encodeURIComponent(gene)}/bins?${query}`);
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

        const queryStr = query.toString();
        const url = `${this.baseUrl}/genes/${encodeURIComponent(gene)}/stats${queryStr ? '?' + queryStr : ''}`;

        return this.fetchJson(url);
    }

    /**
     * Get mean expression per category value for a gene
     */
    async getGeneCategoryMeans(gene: string, column: string): Promise<GeneCategoryMeanItem[]> {
        const data: GeneCategoryMeansResponse = await this.fetchJson(
            `${this.baseUrl}/genes/${encodeURIComponent(gene)}/category/${encodeURIComponent(column)}/means`
        );
        return data.items;
    }

    async getSomaObsColumns(): Promise<string[]> {
        const data: SomaObsColumnsResponse = await this.fetchJson(`${this.baseUrl}/soma/obs/columns`);
        return data.columns;
    }

    async getSomaObsColumnValues(column: string): Promise<string[]> {
        const data: SomaObsValuesResponse = await this.fetchJson(
            `${this.baseUrl}/soma/obs/${encodeURIComponent(column)}/values`
        );
        return data.values;
    }

    async createSomaDeJob(params: SomaDeJobCreateRequest): Promise<SomaDeJobCreateResponse> {
        return this.fetchJson(`${this.baseUrl}/soma/de/jobs`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(params),
        });
    }

    async getSomaDeJob(jobId: string): Promise<SomaDeJobStatusResponse> {
        return this.fetchJson(`${this.baseUrl}/soma/de/jobs/${encodeURIComponent(jobId)}`);
    }

    async getSomaDeJobResult(
        jobId: string,
        params: { offset?: number; limit?: number; order_by?: SomaDeResultOrderBy } = {}
    ): Promise<SomaDeJobResultResponse> {
        const query = new URLSearchParams();
        if (params.offset !== undefined) query.set('offset', params.offset.toString());
        if (params.limit !== undefined) query.set('limit', params.limit.toString());
        if (params.order_by !== undefined) query.set('order_by', params.order_by);
        const suffix = query.toString() ? `?${query}` : '';
        return this.fetchJson(`${this.baseUrl}/soma/de/jobs/${encodeURIComponent(jobId)}/result${suffix}`);
    }

    async cancelSomaDeJob(jobId: string): Promise<void> {
        const response = await fetch(`${this.baseUrl}/soma/de/jobs/${encodeURIComponent(jobId)}`, {
            method: 'DELETE',
        });
        if (!response.ok) {
            const message = await this.readError(response);
            throw new Error(`HTTP ${response.status}: ${message}`);
        }
    }
}
