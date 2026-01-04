// Category Filter Component

import { CategoryInfo } from '../api/client';

export class CategoryFilter {
    private container: HTMLElement;
    private categories: Record<string, CategoryInfo>;
    private selected: Set<string> = new Set();
    private filterCallback: ((categories: string[]) => void) | null = null;

    constructor(container: HTMLElement, categories: Record<string, CategoryInfo>) {
        this.container = container;
        this.categories = categories;
        // Initialize selected with all category values (checkboxes start checked)
        Object.values(categories).forEach(catInfo => {
            catInfo.values.forEach(v => this.selected.add(v));
        });
        this.render();
    }

    private render(): void {
        const categoryNames = Object.keys(this.categories);

        if (categoryNames.length === 0) {
            this.container.innerHTML = '<p class="placeholder">No categories available</p>';
            return;
        }

        // For each category column
        this.container.innerHTML = categoryNames
            .map(colName => {
                const catInfo = this.categories[colName];
                return `
                    <div class="category-group">
                        <h3 class="category-group-title">${colName}</h3>
                        <div class="category-list">
                            ${catInfo.values.map((value, idx) => `
                                <label class="category-item">
                                    <input type="checkbox"
                                           data-column="${colName}"
                                           data-value="${value}"
                                           checked>
                                    <span class="category-color" style="background-color: ${this.getColor(idx)}"></span>
                                    <span class="category-label">${value}</span>
                                </label>
                            `).join('')}
                        </div>
                        <div class="category-actions">
                            <button class="btn-link" data-action="select-all" data-column="${colName}">All</button>
                            <button class="btn-link" data-action="select-none" data-column="${colName}">None</button>
                        </div>
                    </div>
                `;
            })
            .join('');

        this.setupEvents();
    }

    private setupEvents(): void {
        // Checkbox change events
        this.container.querySelectorAll('input[type="checkbox"]').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                const target = e.target as HTMLInputElement;
                const value = target.dataset.value!;

                if (target.checked) {
                    this.selected.add(value);
                } else {
                    this.selected.delete(value);
                }

                this.notifyChange();
            });
        });

        // Select all/none buttons
        this.container.querySelectorAll('button[data-action]').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const target = e.target as HTMLButtonElement;
                const action = target.dataset.action;
                const column = target.dataset.column!;

                const checkboxes = this.container.querySelectorAll(
                    `input[data-column="${column}"]`
                ) as NodeListOf<HTMLInputElement>;

                checkboxes.forEach(cb => {
                    cb.checked = action === 'select-all';
                    const value = cb.dataset.value!;
                    if (action === 'select-all') {
                        this.selected.add(value);
                    } else {
                        this.selected.delete(value);
                    }
                });

                this.notifyChange();
            });
        });
    }

    private getColor(index: number): string {
        // Same colors as Go backend categorical colormap
        const colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
            '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
            '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
        ];
        return colors[index % colors.length];
    }

    private notifyChange(): void {
        if (this.filterCallback) {
            this.filterCallback(Array.from(this.selected));
        }
    }

    /**
     * Register callback for filter changes
     */
    onFilterChange(callback: (categories: string[]) => void): void {
        this.filterCallback = callback;
    }

    /**
     * Get selected categories
     */
    getSelectedCategories(): string[] {
        return Array.from(this.selected);
    }

    /**
     * Select all categories
     */
    selectAll(): void {
        Object.values(this.categories).forEach(catInfo => {
            catInfo.values.forEach(v => this.selected.add(v));
        });
        this.render();
        this.notifyChange();
    }

    /**
     * Clear all selections
     */
    clearAll(): void {
        this.selected.clear();
        this.render();
        this.notifyChange();
    }
}
