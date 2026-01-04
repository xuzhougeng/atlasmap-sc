// Category Legend Component
// Displays the color legend for category coloring

import { ApiClient, CategoryLegendItem } from '../api/client';

export class CategoryLegend {
    private container: HTMLElement;
    private api: ApiClient;
    private items: CategoryLegendItem[] = [];
    private isCollapsed: boolean = false;
    private currentColumn: string = '';

    constructor(container: HTMLElement, api: ApiClient) {
        this.container = container;
        this.api = api;
    }

    /**
     * Load and display legend for a category column
     */
    async loadLegend(column: string): Promise<void> {
        this.currentColumn = column;

        try {
            this.items = await this.api.getCategoryLegend(column);
            this.render();
        } catch (error) {
            console.error('Failed to load legend:', error);
            this.container.innerHTML = '<p class="error">Failed to load legend</p>';
        }
    }

    private render(): void {
        if (this.items.length === 0) {
            this.container.innerHTML = '';
            return;
        }

        this.container.innerHTML = `
            <div class="category-legend">
                <div class="legend-header" id="legend-header">
                    <span class="legend-title">${this.currentColumn}</span>
                    <span class="legend-toggle">${this.isCollapsed ? '+' : '-'}</span>
                </div>
                <div class="legend-items ${this.isCollapsed ? 'collapsed' : ''}">
                    ${this.items.map(item => `
                        <div class="legend-item" data-value="${item.value}" data-index="${item.index}">
                            <span class="legend-color" style="background-color: ${item.color}"></span>
                            <span class="legend-label" title="${item.value}">${item.value}</span>
                            <span class="legend-count">${item.cell_count.toLocaleString()}</span>
                        </div>
                    `).join('')}
                </div>
            </div>
        `;

        this.setupEvents();
    }

    private setupEvents(): void {
        // Toggle collapse/expand
        const header = this.container.querySelector('#legend-header');
        if (header) {
            header.addEventListener('click', () => {
                this.isCollapsed = !this.isCollapsed;
                this.render();
            });
        }

        // Hover effect on legend items
        this.container.querySelectorAll('.legend-item').forEach(item => {
            item.addEventListener('mouseenter', () => {
                item.classList.add('highlighted');
            });
            item.addEventListener('mouseleave', () => {
                item.classList.remove('highlighted');
            });
        });
    }

    /**
     * Clear the legend
     */
    clear(): void {
        this.items = [];
        this.currentColumn = '';
        this.container.innerHTML = '';
    }

    /**
     * Show the legend
     */
    show(): void {
        this.container.style.display = 'block';
    }

    /**
     * Hide the legend
     */
    hide(): void {
        this.container.style.display = 'none';
    }

    /**
     * Get current column name
     */
    getCurrentColumn(): string {
        return this.currentColumn;
    }
}
