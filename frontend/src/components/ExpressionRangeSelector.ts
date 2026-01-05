// Expression Range Selector Component
// Controls min/max range for expression colormap scaling

export interface ExpressionRange {
    min: number;
    max: number;
}

export interface ExpressionRangeSelectorConfig {
    onRangeChange: (range: ExpressionRange | null) => void;
}

export class ExpressionRangeSelector {
    private container: HTMLElement;
    private config: ExpressionRangeSelectorConfig;
    private auto: boolean = true;
    private minValue: number | null = null;
    private maxValue: number | null = null;

    constructor(container: HTMLElement, config: ExpressionRangeSelectorConfig) {
        this.container = container;
        this.config = config;
        this.render();
    }

    private render(): void {
        this.container.innerHTML = `
            <div class="expression-range-selector">
                <div class="expression-range-header">
                    <label>Range:</label>
                    <label class="expression-range-auto">
                        <input type="checkbox" id="expression-range-auto" ${this.auto ? 'checked' : ''}>
                        Auto
                    </label>
                </div>
                <div class="expression-range-inputs">
                    <input type="number"
                           id="expression-range-min"
                           class="expression-range-input"
                           placeholder="min"
                           step="any"
                           ${this.auto ? 'disabled' : ''}
                           value="${this.minValue ?? ''}">
                    <span class="expression-range-sep">–</span>
                    <input type="number"
                           id="expression-range-max"
                           class="expression-range-input"
                           placeholder="max"
                           step="any"
                           ${this.auto ? 'disabled' : ''}
                           value="${this.maxValue ?? ''}">
                </div>
            </div>
        `;

        this.setupEvents();
    }

    private setupEvents(): void {
        const autoEl = this.container.querySelector('#expression-range-auto') as HTMLInputElement | null;
        const minEl = this.container.querySelector('#expression-range-min') as HTMLInputElement | null;
        const maxEl = this.container.querySelector('#expression-range-max') as HTMLInputElement | null;
        if (!autoEl || !minEl || !maxEl) return;

        autoEl.addEventListener('change', () => {
            this.auto = autoEl.checked;
            if (this.auto) {
                minEl.disabled = true;
                maxEl.disabled = true;
                this.config.onRangeChange(null);
                return;
            }

            minEl.disabled = false;
            maxEl.disabled = false;
            if (this.minValue === null) {
                this.minValue = 0;
                minEl.value = '0';
            }
            if (this.maxValue === null) {
                this.maxValue = 1;
                maxEl.value = '1';
            }
            const range = this.getRange();
            if (range) {
                this.config.onRangeChange(range);
            }
        });

        const onInputChange = () => {
            this.minValue = this.parseNumber(minEl.value);
            this.maxValue = this.parseNumber(maxEl.value);
            const range = this.getRange();
            if (range) {
                this.config.onRangeChange(range);
            }
        };

        minEl.addEventListener('change', onInputChange);
        maxEl.addEventListener('change', onInputChange);

        [minEl, maxEl].forEach(el => {
            el.addEventListener('keydown', (e: KeyboardEvent) => {
                if (e.key === 'Enter') {
                    (e.target as HTMLInputElement).blur();
                }
            });
        });
    }

    private parseNumber(value: string): number | null {
        const trimmed = value.trim();
        if (!trimmed) return null;
        const parsed = Number(trimmed);
        if (!Number.isFinite(parsed)) return null;
        return parsed;
    }

    getRange(): ExpressionRange | null {
        if (this.auto) return null;
        if (this.minValue === null || this.maxValue === null) return null;
        const min = Math.min(this.minValue, this.maxValue);
        const max = Math.max(this.minValue, this.maxValue);
        return { min, max };
    }
}
