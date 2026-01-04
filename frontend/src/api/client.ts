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

export class ApiClient {
    private baseUrl: string;

    constructor(baseUrl: string) {
        this.baseUrl = baseUrl;
    }

    async getMetadata(): Promise<Metadata> {
        const response = await fetch(`${this.baseUrl}/metadata`);
        if (!response.ok) {
            throw new Error(`Failed to fetch metadata: ${response.statusText}`);
        }
        return response.json();
    }

    async getGenes(): Promise<{ genes: string[]; total: number }> {
        const response = await fetch(`${this.baseUrl}/genes`);
        if (!response.ok) {
            throw new Error(`Failed to fetch genes: ${response.statusText}`);
        }
        return response.json();
    }

    async getGeneInfo(gene: string): Promise<GeneInfo> {
        const response = await fetch(`${this.baseUrl}/genes/${encodeURIComponent(gene)}`);
        if (!response.ok) {
            throw new Error(`Failed to fetch gene info: ${response.statusText}`);
        }
        return response.json();
    }

    async getCategories(): Promise<Record<string, CategoryInfo>> {
        const response = await fetch(`${this.baseUrl}/categories`);
        if (!response.ok) {
            throw new Error(`Failed to fetch categories: ${response.statusText}`);
        }
        return response.json();
    }

    async getStats(): Promise<Record<string, any>> {
        const response = await fetch(`${this.baseUrl}/stats`);
        if (!response.ok) {
            throw new Error(`Failed to fetch stats: ${response.statusText}`);
        }
        return response.json();
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
        const response = await fetch(
            `${this.baseUrl}/categories/${encodeURIComponent(column)}/colors`
        );
        if (!response.ok) {
            throw new Error(`Failed to fetch category colors: ${response.statusText}`);
        }
        return response.json();
    }

    /**
     * Get legend data for a category column
     */
    async getCategoryLegend(column: string): Promise<CategoryLegendItem[]> {
        const response = await fetch(
            `${this.baseUrl}/categories/${encodeURIComponent(column)}/legend`
        );
        if (!response.ok) {
            throw new Error(`Failed to fetch category legend: ${response.statusText}`);
        }
        return response.json();
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

        const response = await fetch(
            `${this.baseUrl}/genes/${encodeURIComponent(gene)}/bins?${query}`
        );
        if (!response.ok) {
            throw new Error(`Failed to query bins: ${response.statusText}`);
        }
        return response.json();
    }

    /**
     * Get gene expression statistics
     */
    async getGeneStats(gene: string): Promise<GeneStats> {
        const response = await fetch(
            `${this.baseUrl}/genes/${encodeURIComponent(gene)}/stats`
        );
        if (!response.ok) {
            throw new Error(`Failed to get gene stats: ${response.statusText}`);
        }
        return response.json();
    }

    /**
     * Get mean expression per category value for a gene
     */
    async getGeneCategoryMeans(gene: string, column: string): Promise<GeneCategoryMeanItem[]> {
        const response = await fetch(
            `${this.baseUrl}/genes/${encodeURIComponent(gene)}/category/${encodeURIComponent(column)}/means`
        );
        if (!response.ok) {
            throw new Error(`Failed to get gene category means: ${response.statusText}`);
        }
        const data: GeneCategoryMeansResponse = await response.json();
        return data.items;
    }
}
