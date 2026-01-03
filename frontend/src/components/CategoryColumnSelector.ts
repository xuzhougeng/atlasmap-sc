// Category Column Selector Component
// Dropdown for selecting the category column to color by

export interface CategoryColumnSelectorConfig {
    onCategoryChange: (column: string) => void;
}

export class CategoryColumnSelector {
    private container: HTMLElement;
    private config: CategoryColumnSelectorConfig;
    private availableCategories: string[] = [];
    private selectedCategory: string = '';

    constructor(container: HTMLElement, config: CategoryColumnSelectorConfig) {
        this.container = container;
        this.config = config;
        this.render();
    }

    private render(): void {
        this.container.innerHTML = `
            <div class="category-column-selector">
                <label>Color by:</label>
                <select id="category-column-select">
                    ${this.availableCategories.map(cat =>
                        `<option value="${cat}" ${cat === this.selectedCategory ? 'selected' : ''}>
                            ${cat}
                        </option>`
                    ).join('')}
                </select>
            </div>
        `;

        this.setupEvents();
    }

    private setupEvents(): void {
        const select = this.container.querySelector('#category-column-select') as HTMLSelectElement;
        if (select) {
            select.addEventListener('change', () => {
                this.selectedCategory = select.value;
                this.config.onCategoryChange(this.selectedCategory);
            });
        }
    }

    setCategories(categories: string[]): void {
        this.availableCategories = categories;
        if (categories.length > 0 && !categories.includes(this.selectedCategory)) {
            this.selectedCategory = categories[0];
        }
        this.render();
    }

    getSelectedCategory(): string {
        return this.selectedCategory;
    }

    setSelectedCategory(category: string): void {
        if (this.availableCategories.includes(category)) {
            this.selectedCategory = category;
            const select = this.container.querySelector('#category-column-select') as HTMLSelectElement;
            if (select) {
                select.value = category;
            }
        }
    }
}
