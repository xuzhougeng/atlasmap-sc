// Tab Panel Component
// Manages tab switching between Category and Expression views

export type TabType = 'category' | 'expression';

export interface TabPanelConfig {
    onTabChange: (tab: TabType) => void;
}

export class TabPanel {
    private container: HTMLElement;
    private config: TabPanelConfig;
    private currentTab: TabType = 'category';

    // Content containers for each tab
    private categoryContent!: HTMLElement;
    private expressionContent!: HTMLElement;

    constructor(container: HTMLElement, config: TabPanelConfig) {
        this.container = container;
        this.config = config;
        this.render();
    }

    private render(): void {
        this.container.innerHTML = `
            <div class="tab-panel">
                <div class="tab-header">
                    <button class="tab-btn ${this.currentTab === 'category' ? 'active' : ''}"
                            data-tab="category">
                        Category
                    </button>
                    <button class="tab-btn ${this.currentTab === 'expression' ? 'active' : ''}"
                            data-tab="expression">
                        Expression
                    </button>
                </div>
                <div class="tab-content">
                    <div id="tab-category" class="tab-pane ${this.currentTab === 'category' ? 'active' : ''}">
                    </div>
                    <div id="tab-expression" class="tab-pane ${this.currentTab === 'expression' ? 'active' : ''}">
                    </div>
                </div>
            </div>
        `;

        this.categoryContent = this.container.querySelector('#tab-category')!;
        this.expressionContent = this.container.querySelector('#tab-expression')!;

        this.setupEvents();
    }

    private setupEvents(): void {
        const tabButtons = this.container.querySelectorAll('.tab-btn');
        tabButtons.forEach(btn => {
            btn.addEventListener('click', (e) => {
                const tab = (e.target as HTMLElement).dataset.tab as TabType;
                this.switchTab(tab);
            });
        });
    }

    switchTab(tab: TabType): void {
        if (this.currentTab === tab) return;

        this.currentTab = tab;

        // Update button styles
        const tabButtons = this.container.querySelectorAll('.tab-btn');
        tabButtons.forEach(btn => {
            btn.classList.toggle('active', (btn as HTMLElement).dataset.tab === tab);
        });

        // Update pane visibility
        this.categoryContent.classList.toggle('active', tab === 'category');
        this.expressionContent.classList.toggle('active', tab === 'expression');

        // Notify parent
        this.config.onTabChange(tab);
    }

    getCurrentTab(): TabType {
        return this.currentTab;
    }

    getCategoryContainer(): HTMLElement {
        return this.categoryContent;
    }

    getExpressionContainer(): HTMLElement {
        return this.expressionContent;
    }
}
