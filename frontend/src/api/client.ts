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
}
