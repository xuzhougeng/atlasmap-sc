// Gene Selector Component

import { ApiClient } from '../api/client';

export class GeneSelector {
    private container: HTMLElement;
    private input: HTMLInputElement;
    private dropdown: HTMLElement;
    private api: ApiClient;
    private genes: string[] = [];
    private selectedGene: string | null = null;
    private selectCallback: ((gene: string) => void) | null = null;
    private debounceTimer: number | null = null;

    constructor(container: HTMLElement, api: ApiClient) {
        this.container = container;
        this.api = api;

        // Get existing elements
        this.input = container.querySelector('#gene-input') as HTMLInputElement;
        this.dropdown = container.querySelector('#gene-dropdown') as HTMLElement;

        if (!this.input || !this.dropdown) {
            throw new Error('Gene selector elements not found');
        }

        this.setupEvents();
        this.loadGenes();
    }

    private async loadGenes(): Promise<void> {
        try {
            const response = await this.api.getGenes();
            this.genes = response.genes;
            console.log(`Loaded ${this.genes.length} genes`);
        } catch (error) {
            console.error('Failed to load genes:', error);
        }
    }

    private setupEvents(): void {
        // Input event for search
        this.input.addEventListener('input', () => {
            this.debounce(() => this.search(this.input.value), 200);
        });

        // Focus event to show dropdown
        this.input.addEventListener('focus', () => {
            if (this.input.value) {
                this.search(this.input.value);
            }
        });

        // Click outside to close dropdown
        document.addEventListener('click', (e) => {
            if (!this.container.contains(e.target as Node)) {
                this.hideDropdown();
            }
        });

        // Keyboard navigation
        this.input.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                this.hideDropdown();
            } else if (e.key === 'Enter') {
                const firstItem = this.dropdown.querySelector('.gene-item') as HTMLElement;
                if (firstItem) {
                    firstItem.click();
                }
            }
        });
    }

    private debounce(fn: () => void, delay: number): void {
        if (this.debounceTimer) {
            clearTimeout(this.debounceTimer);
        }
        this.debounceTimer = window.setTimeout(fn, delay);
    }

    private search(query: string): void {
        if (!query || query.length < 1) {
            this.hideDropdown();
            return;
        }

        const queryLower = query.toLowerCase();
        const matches = this.genes
            .filter(g => g.toLowerCase().includes(queryLower))
            .slice(0, 20);

        if (matches.length === 0) {
            this.hideDropdown();
            return;
        }

        this.showDropdown(matches);
    }

    private showDropdown(genes: string[]): void {
        this.dropdown.innerHTML = genes
            .map(gene => `
                <div class="gene-item" data-gene="${gene}">
                    ${this.highlightMatch(gene, this.input.value)}
                </div>
            `)
            .join('');

        // Add click handlers
        this.dropdown.querySelectorAll('.gene-item').forEach(item => {
            item.addEventListener('click', () => {
                const gene = (item as HTMLElement).dataset.gene!;
                this.selectGene(gene);
            });
        });

        this.dropdown.classList.add('visible');
    }

    private hideDropdown(): void {
        this.dropdown.classList.remove('visible');
    }

    private highlightMatch(gene: string, query: string): string {
        const idx = gene.toLowerCase().indexOf(query.toLowerCase());
        if (idx === -1) return gene;

        return (
            gene.slice(0, idx) +
            `<strong>${gene.slice(idx, idx + query.length)}</strong>` +
            gene.slice(idx + query.length)
        );
    }

    private selectGene(gene: string): void {
        this.selectedGene = gene;
        this.input.value = gene;
        this.hideDropdown();

        if (this.selectCallback) {
            this.selectCallback(gene);
        }
    }

    /**
     * Register callback for gene selection
     */
    onGeneSelect(callback: (gene: string) => void): void {
        this.selectCallback = callback;
    }

    /**
     * Get currently selected gene
     */
    getSelectedGene(): string | null {
        return this.selectedGene;
    }

    /**
     * Clear selection
     */
    clear(): void {
        this.selectedGene = null;
        this.input.value = '';
    }
}
