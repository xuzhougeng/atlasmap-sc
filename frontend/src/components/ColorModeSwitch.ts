// Color Mode Switch Component
// Allows switching between expression and category coloring modes

export type ColorMode = 'expression' | 'category';

export interface ColorModeSwitchConfig {
    onModeChange: (mode: ColorMode) => void;
    onCategoryChange?: (column: string) => void;
}

export class ColorModeSwitch {
    private container: HTMLElement;
    private config: ColorModeSwitchConfig;
    private currentMode: ColorMode = 'category';
    private availableCategories: string[] = [];
    private selectedCategory: string = 'cell_type';

    constructor(container: HTMLElement, config: ColorModeSwitchConfig) {
        this.container = container;
        this.config = config;
        this.render();
    }

    private render(): void {
        this.container.innerHTML = `
            <div class="color-mode-switch">
                <div class="mode-buttons">
                    <button id="btn-mode-category"
                            class="mode-btn ${this.currentMode === 'category' ? 'active' : ''}">
                        Category
                    </button>
                    <button id="btn-mode-expression"
                            class="mode-btn ${this.currentMode === 'expression' ? 'active' : ''}">
                        Expression
                    </button>
                </div>

                <div id="category-selector"
                     class="category-selector ${this.currentMode === 'category' ? '' : 'hidden'}">
                    <label>Color by:</label>
                    <select id="category-column-select">
                        ${this.availableCategories.map(cat =>
                            `<option value="${cat}" ${cat === this.selectedCategory ? 'selected' : ''}>
                                ${cat}
                            </option>`
                        ).join('')}
                    </select>
                </div>
            </div>
        `;

        this.setupEvents();
    }

    private setupEvents(): void {
        // Mode switch buttons
        const exprBtn = this.container.querySelector('#btn-mode-expression');
        const catBtn = this.container.querySelector('#btn-mode-category');

        if (exprBtn) {
            exprBtn.addEventListener('click', () => this.setMode('expression'));
        }
        if (catBtn) {
            catBtn.addEventListener('click', () => this.setMode('category'));
        }

        // Category column select
        const select = this.container.querySelector('#category-column-select') as HTMLSelectElement;
        if (select) {
            select.addEventListener('change', () => {
                this.selectedCategory = select.value;
                if (this.config.onCategoryChange) {
                    this.config.onCategoryChange(this.selectedCategory);
                }
            });
        }
    }

    private setMode(mode: ColorMode): void {
        if (this.currentMode === mode) return;

        this.currentMode = mode;
        this.config.onModeChange(mode);

        // Update UI
        const exprBtn = this.container.querySelector('#btn-mode-expression');
        const catBtn = this.container.querySelector('#btn-mode-category');
        const catSelector = this.container.querySelector('#category-selector');

        if (exprBtn) {
            exprBtn.classList.toggle('active', mode === 'expression');
        }
        if (catBtn) {
            catBtn.classList.toggle('active', mode === 'category');
        }
        if (catSelector) {
            catSelector.classList.toggle('hidden', mode !== 'category');
        }
    }

    /**
     * Set available category columns
     */
    setCategories(categories: string[]): void {
        this.availableCategories = categories;
        if (categories.length > 0 && !categories.includes(this.selectedCategory)) {
            this.selectedCategory = categories[0];
        }
        this.render();
    }

    /**
     * Get current mode
     */
    getMode(): ColorMode {
        return this.currentMode;
    }

    /**
     * Get selected category column
     */
    getSelectedCategory(): string {
        return this.selectedCategory;
    }

    /**
     * Programmatically set mode (e.g., when gene is selected)
     */
    switchToExpression(): void {
        this.setMode('expression');
    }

    /**
     * Programmatically set mode back to category
     */
    switchToCategory(): void {
        this.setMode('category');
    }
}
