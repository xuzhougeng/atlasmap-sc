// Category Filter Component

import { CategoryInfo } from '../api/client';

export class CategoryFilter {
    private container: HTMLElement;
    private categories: Record<string, CategoryInfo>;
    private selectedByColumn: Map<string, Set<string>> = new Map();
    private filterCallback: ((column: string, categories: string[] | null) => void) | null = null;
    private activeColumn: string | null = null;

    constructor(
        container: HTMLElement,
        categories: Record<string, CategoryInfo>,
        activeColumn?: string
    ) {
        this.container = container;
        this.categories = categories;
        this.activeColumn =
            activeColumn && Object.prototype.hasOwnProperty.call(categories, activeColumn)
                ? activeColumn
                : null;
        // Initialize selected sets with all category values (checkboxes start checked)
        for (const [column, catInfo] of Object.entries(categories)) {
            this.selectedByColumn.set(column, new Set(catInfo.values));
        }
        this.render();
    }

    private render(): void {
        const allColumns = Object.keys(this.categories);
        const columnsToRender = this.activeColumn ? [this.activeColumn] : allColumns;

        if (columnsToRender.length === 0) {
            this.container.innerHTML = '<p class="placeholder">No categories available</p>';
            return;
        }

        const validColumns = columnsToRender.filter(colName =>
            Object.prototype.hasOwnProperty.call(this.categories, colName)
        );

        if (validColumns.length === 0) {
            this.container.innerHTML = '<p class="placeholder">Select a category to filter</p>';
            return;
        }

        // For each (visible) category column
        this.container.innerHTML = validColumns
            .map(colName => {
                const catInfo = this.categories[colName];
                const selected = this.selectedByColumn.get(colName) ?? new Set<string>();
                return `
                    <div class="category-group">
                        <h3 class="category-group-title">${colName}</h3>
                        <div class="category-list">
                            ${catInfo.values.map((value, idx) => `
                                <label class="category-item">
                                    <input type="checkbox"
                                           data-column="${colName}"
                                           data-value="${value}"
                                           ${selected.has(value) ? 'checked' : ''}>
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

    setActiveColumn(column: string): void {
        const next =
            column && Object.prototype.hasOwnProperty.call(this.categories, column) ? column : null;
        if (this.activeColumn === next) return;
        this.activeColumn = next;
        this.render();
    }

    private setupEvents(): void {
        // Checkbox change events
        this.container.querySelectorAll('input[type="checkbox"]').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                const target = e.target as HTMLInputElement;
                const column = target.dataset.column!;
                const value = target.dataset.value!;
                const selected = this.selectedByColumn.get(column);
                if (!selected) return;

                if (target.checked) {
                    selected.add(value);
                } else {
                    selected.delete(value);
                }

                this.notifyChange(column);
            });
        });

        // Select all/none buttons
        this.container.querySelectorAll('button[data-action]').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const target = e.target as HTMLButtonElement;
                const action = target.dataset.action;
                const column = target.dataset.column!;
                const selected = this.selectedByColumn.get(column);
                if (!selected) return;

                const checkboxes = this.container.querySelectorAll(
                    `input[data-column="${column}"]`
                ) as NodeListOf<HTMLInputElement>;

                checkboxes.forEach(cb => {
                    cb.checked = action === 'select-all';
                    const value = cb.dataset.value!;
                    if (action === 'select-all') {
                        selected.add(value);
                    } else {
                        selected.delete(value);
                    }
                });

                this.notifyChange(column);
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

    private notifyChange(column: string): void {
        if (!this.filterCallback) return;
        this.filterCallback(column, this.getFilterForColumn(column));
    }

    /**
     * Register callback for filter changes
     */
    onFilterChange(callback: (column: string, categories: string[] | null) => void): void {
        this.filterCallback = callback;
    }

    /**
     * Get column filter.
     *
     * Returns null when all values are selected (no filter).
     * Returns an array (possibly empty) when a filter is active.
     */
    getFilterForColumn(column: string): string[] | null {
        const catInfo = this.categories[column];
        const selected = this.selectedByColumn.get(column);
        if (!catInfo || !selected) return null;
        if (selected.size === catInfo.values.length) return null;
        return Array.from(selected);
    }

    /**
     * Select all categories
     */
    selectAll(): void {
        for (const [column, catInfo] of Object.entries(this.categories)) {
            this.selectedByColumn.set(column, new Set(catInfo.values));
        }
        this.render();
        for (const column of Object.keys(this.categories)) {
            this.notifyChange(column);
        }
    }

    /**
     * Clear all selections
     */
    clearAll(): void {
        for (const column of Object.keys(this.categories)) {
            this.selectedByColumn.set(column, new Set());
        }
        this.render();
        for (const column of Object.keys(this.categories)) {
            this.notifyChange(column);
        }
    }
}
