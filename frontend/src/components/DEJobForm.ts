// Differential Expression Job Form Component
// Allows submitting DE analysis jobs

import { ApiClient, SomaDeTest, SomaDeJobCreateRequest } from '../api/client';

export interface DEJobFormConfig {
    /** Called when job is successfully submitted */
    onSubmit?: (jobId: string) => void;
    /** Called when submission fails */
    onError?: (error: Error) => void;
}

export class DEJobForm {
    private container: HTMLElement;
    private api: ApiClient;
    private config: DEJobFormConfig;

    private allGroupValues: string[] = [];
    private isSubmitting: boolean = false;

    constructor(
        container: HTMLElement,
        api: ApiClient,
        config: DEJobFormConfig = {}
    ) {
        this.container = container;
        this.api = api;
        this.config = config;
        this.render();
        this.loadColumns();
    }

    private render(): void {
        this.container.innerHTML = `
            <div class="de-job-form">
                <div class="form-grid">
                    <div class="form-row">
                        <label>Test Methods</label>
                        <div class="form-field">
                            <label class="checkbox-label">
                                <input type="checkbox" id="de-test-ttest" checked /> t-test
                            </label>
                            <label class="checkbox-label">
                                <input type="checkbox" id="de-test-ranksum" checked /> rank-sum
                            </label>
                            <div class="form-help">Welch t-test / Mann-Whitney U (Wilcoxon rank-sum)</div>
                        </div>
                    </div>

                    <div class="form-row">
                        <label for="de-groupby">Group By</label>
                        <div class="form-field">
                            <select id="de-groupby" class="select" disabled>
                                <option value="">Loading...</option>
                            </select>
                            <div class="form-help">Column from obs metadata</div>
                        </div>
                    </div>

                    <div class="form-row">
                        <label>Group Selection</label>
                        <div class="form-field">
                            <div class="group-selection-controls">
                                <input id="de-value-filter" class="input" placeholder="Filter values (optional)" />
                                <label class="checkbox-label">
                                    <input type="checkbox" id="de-one-vs-rest" />
                                    One-vs-rest (leave group2 empty)
                                </label>
                            </div>
                            <div class="multi-select">
                                <div class="multi-select-box">
                                    <div class="multi-select-title">Group 1 (at least 1)</div>
                                    <div class="multi-select-actions">
                                        <button class="btn-link" id="de-g1-all" type="button">Select All</button>
                                        <button class="btn-link" id="de-g1-clear" type="button">Clear</button>
                                    </div>
                                    <select id="de-group1" class="select" multiple size="10" disabled></select>
                                    <div class="chip-list" id="de-group1-chips"></div>
                                </div>
                                <div class="multi-select-box">
                                    <div class="multi-select-title">Group 2 (optional)</div>
                                    <div class="multi-select-actions">
                                        <button class="btn-link" id="de-g2-all" type="button">Select All</button>
                                        <button class="btn-link" id="de-g2-clear" type="button">Clear</button>
                                    </div>
                                    <select id="de-group2" class="select" multiple size="10" disabled></select>
                                    <div class="chip-list" id="de-group2-chips"></div>
                                </div>
                            </div>
                            <div class="form-help">Leave group2 empty for one-vs-rest comparison (group1 vs all others).</div>
                        </div>
                    </div>

                    <div class="form-row">
                        <label>Advanced Options</label>
                        <div class="form-field">
                            <div class="param-grid">
                                <div class="param-item">
                                    <label class="param-label">max_cells_per_group</label>
                                    <input id="de-max-cells" class="input" type="number" min="1" max="20000" value="2000" />
                                </div>
                                <div class="param-item">
                                    <label class="param-label">seed</label>
                                    <input id="de-seed" class="input" type="number" value="0" />
                                </div>
                                <div class="param-item">
                                    <label class="param-label">limit</label>
                                    <input id="de-limit" class="input" type="number" min="1" max="500" value="50" />
                                </div>
                            </div>
                        </div>
                    </div>

                    <div class="form-row">
                        <label></label>
                        <div class="form-field form-actions">
                            <button id="de-btn-submit" class="btn-primary" type="button" disabled>Submit DE Job</button>
                            <span id="de-submit-hint" class="form-help"></span>
                        </div>
                    </div>

                    <div id="de-form-error"></div>
                </div>
            </div>
        `;

        this.setupEvents();
    }

    private setupEvents(): void {
        // Test method checkboxes
        this.container.querySelectorAll('#de-test-ttest, #de-test-ranksum').forEach(el => {
            el.addEventListener('change', () => this.updateSubmitState());
        });

        // Group by select
        const groupbyEl = this.getElement<HTMLSelectElement>('#de-groupby');
        groupbyEl.addEventListener('change', () => this.handleGroupbyChange());

        // Filter input
        const filterEl = this.getElement<HTMLInputElement>('#de-value-filter');
        filterEl.addEventListener('input', () => {
            this.applyFilter();
            this.updateChips();
        });

        // One-vs-rest checkbox
        const oneVsRestEl = this.getElement<HTMLInputElement>('#de-one-vs-rest');
        oneVsRestEl.addEventListener('change', () => {
            const group2El = this.getElement<HTMLSelectElement>('#de-group2');
            group2El.disabled = oneVsRestEl.checked || !groupbyEl.value;
            this.updateSubmitState();
        });

        // Group selects
        const group1El = this.getElement<HTMLSelectElement>('#de-group1');
        const group2El = this.getElement<HTMLSelectElement>('#de-group2');
        group1El.addEventListener('change', () => {
            this.updateChips();
            this.updateSubmitState();
        });
        group2El.addEventListener('change', () => {
            this.updateChips();
            this.updateSubmitState();
        });

        // Select All / Clear buttons
        this.getElement('#de-g1-all').addEventListener('click', () => {
            this.selectAll(group1El);
            this.updateChips();
            this.updateSubmitState();
        });
        this.getElement('#de-g1-clear').addEventListener('click', () => {
            this.clearAll(group1El);
            this.updateChips();
            this.updateSubmitState();
        });
        this.getElement('#de-g2-all').addEventListener('click', () => {
            this.selectAll(group2El);
            this.updateChips();
            this.updateSubmitState();
        });
        this.getElement('#de-g2-clear').addEventListener('click', () => {
            this.clearAll(group2El);
            this.updateChips();
            this.updateSubmitState();
        });

        // Submit button
        this.getElement('#de-btn-submit').addEventListener('click', () => this.handleSubmit());
    }

    private getElement<T extends HTMLElement>(selector: string): T {
        return this.container.querySelector(selector) as T;
    }

    private async loadColumns(): Promise<void> {
        const groupbyEl = this.getElement<HTMLSelectElement>('#de-groupby');

        try {
            const columns = await this.api.getSomaObsColumns();
            groupbyEl.innerHTML = '<option value="">Select...</option>';
            for (const c of columns) {
                const opt = document.createElement('option');
                opt.value = c;
                opt.textContent = c;
                groupbyEl.appendChild(opt);
            }
            groupbyEl.disabled = false;
        } catch (e: any) {
            console.error('Failed to load obs columns:', e);
            this.setError(`Failed to load obs columns: ${e?.message || String(e)} (requires backend built with -tags soma)`);
        }

        this.updateSubmitState();
    }

    private async handleGroupbyChange(): Promise<void> {
        this.setError(null);
        this.allGroupValues = [];

        const group1El = this.getElement<HTMLSelectElement>('#de-group1');
        const group2El = this.getElement<HTMLSelectElement>('#de-group2');
        const oneVsRestEl = this.getElement<HTMLInputElement>('#de-one-vs-rest');
        const groupbyEl = this.getElement<HTMLSelectElement>('#de-groupby');

        group1El.innerHTML = '';
        group2El.innerHTML = '';
        group1El.disabled = true;
        group2El.disabled = true;
        this.updateChips();
        this.updateSubmitState();

        const column = groupbyEl.value;
        if (!column) return;

        group1El.disabled = false;
        group2El.disabled = oneVsRestEl.checked;
        group1El.innerHTML = '<option>Loading...</option>';
        group2El.innerHTML = '<option>Loading...</option>';

        try {
            this.allGroupValues = this.uniqSorted(await this.api.getSomaObsColumnValues(column));
        } catch (e: any) {
            console.error('Failed to load group values:', e);
            this.setError(`Failed to load group values: ${e?.message || String(e)}`);
            group1El.innerHTML = '';
            group2El.innerHTML = '';
            group1El.disabled = true;
            group2El.disabled = true;
            this.updateSubmitState();
            return;
        }

        this.applyFilter();
        this.updateChips();
        this.updateSubmitState();
    }

    private applyFilter(): void {
        const filterEl = this.getElement<HTMLInputElement>('#de-value-filter');
        const group1El = this.getElement<HTMLSelectElement>('#de-group1');
        const group2El = this.getElement<HTMLSelectElement>('#de-group2');

        const q = filterEl.value.trim().toLowerCase();
        const filtered = q
            ? this.allGroupValues.filter(v => v.toLowerCase().includes(q))
            : this.allGroupValues;

        const selected1 = new Set(this.getSelectedValues('#de-group1'));
        const selected2 = new Set(this.getSelectedValues('#de-group2'));

        this.renderOptions(group1El, filtered, selected1);
        this.renderOptions(group2El, filtered, selected2);
    }

    private renderOptions(select: HTMLSelectElement, values: string[], selected: Set<string>): void {
        select.innerHTML = '';
        for (const v of values) {
            const opt = document.createElement('option');
            opt.value = v;
            opt.textContent = v;
            opt.selected = selected.has(v);
            select.appendChild(opt);
        }
    }

    private getSelectedValues(selector: string): string[] {
        const select = this.getElement<HTMLSelectElement>(selector);
        return Array.from(select.selectedOptions).map(o => o.value);
    }

    private setSelected(select: HTMLSelectElement, selected: Set<string>): void {
        for (const option of Array.from(select.options)) {
            option.selected = selected.has(option.value);
        }
    }

    private selectAll(select: HTMLSelectElement): void {
        for (const opt of Array.from(select.options)) {
            opt.selected = true;
        }
    }

    private clearAll(select: HTMLSelectElement): void {
        for (const opt of Array.from(select.options)) {
            opt.selected = false;
        }
    }

    private updateChips(): void {
        const selected1 = this.uniqSorted(this.getSelectedValues('#de-group1'));
        const selected2 = this.uniqSorted(this.getSelectedValues('#de-group2'));

        this.renderChips(
            this.getElement('#de-group1-chips'),
            selected1,
            (v) => {
                const set = new Set(this.getSelectedValues('#de-group1'));
                set.delete(v);
                this.setSelected(this.getElement<HTMLSelectElement>('#de-group1'), set);
                this.updateChips();
                this.updateSubmitState();
            }
        );

        this.renderChips(
            this.getElement('#de-group2-chips'),
            selected2,
            (v) => {
                const set = new Set(this.getSelectedValues('#de-group2'));
                set.delete(v);
                this.setSelected(this.getElement<HTMLSelectElement>('#de-group2'), set);
                this.updateChips();
                this.updateSubmitState();
            }
        );
    }

    private renderChips(container: HTMLElement, values: string[], onRemove: (v: string) => void): void {
        container.innerHTML = '';
        for (const v of values) {
            const chip = document.createElement('span');
            chip.className = 'chip';
            chip.innerHTML = `<span>${this.escapeHtml(v)}</span>`;

            const btn = document.createElement('button');
            btn.type = 'button';
            btn.title = 'Remove';
            btn.textContent = '\u00d7';
            btn.addEventListener('click', () => onRemove(v));

            chip.appendChild(btn);
            container.appendChild(chip);
        }
    }

    private updateSubmitState(): void {
        const tests: SomaDeTest[] = [];
        if (this.getElement<HTMLInputElement>('#de-test-ttest').checked) tests.push('ttest');
        if (this.getElement<HTMLInputElement>('#de-test-ranksum').checked) tests.push('ranksum');

        const groupby = this.getElement<HTMLSelectElement>('#de-groupby').value;
        const group1 = this.getSelectedValues('#de-group1');
        const oneVsRest = this.getElement<HTMLInputElement>('#de-one-vs-rest').checked;
        const group2 = oneVsRest ? [] : this.getSelectedValues('#de-group2');

        const isValid = Boolean(groupby) && group1.length > 0 && tests.length > 0;

        this.getElement<HTMLButtonElement>('#de-btn-submit').disabled = !isValid || this.isSubmitting;

        // Update hint
        const hints: string[] = [];
        if (!groupby) hints.push('Select a group by column');
        if (group1.length === 0) hints.push('Select at least one value for group1');
        if (tests.length === 0) hints.push('Select at least one test method');
        if (!oneVsRest && group2.length === 0) hints.push('Group2 empty = one-vs-rest');

        this.getElement('#de-submit-hint').textContent = hints.join('; ');
    }

    private async handleSubmit(): Promise<void> {
        if (this.isSubmitting) return;

        this.isSubmitting = true;
        const submitBtn = this.getElement<HTMLButtonElement>('#de-btn-submit');
        submitBtn.disabled = true;
        submitBtn.textContent = 'Submitting...';
        this.setError(null);

        try {
            const params = this.buildParams();
            const result = await this.api.createSomaDeJob(params);
            this.config.onSubmit?.(result.job_id);
        } catch (e: any) {
            const error = e instanceof Error ? e : new Error(String(e));
            this.setError(`Submission failed: ${error.message}`);
            this.config.onError?.(error);
        } finally {
            this.isSubmitting = false;
            submitBtn.textContent = 'Submit DE Job';
            this.updateSubmitState();
        }
    }

    private buildParams(): SomaDeJobCreateRequest {
        const tests: SomaDeTest[] = [];
        if (this.getElement<HTMLInputElement>('#de-test-ttest').checked) tests.push('ttest');
        if (this.getElement<HTMLInputElement>('#de-test-ranksum').checked) tests.push('ranksum');

        const groupby = this.getElement<HTMLSelectElement>('#de-groupby').value;
        const group1 = this.getSelectedValues('#de-group1');
        const oneVsRest = this.getElement<HTMLInputElement>('#de-one-vs-rest').checked;
        const group2 = oneVsRest ? [] : this.getSelectedValues('#de-group2');

        const maxCells = parseInt(this.getElement<HTMLInputElement>('#de-max-cells').value, 10);
        const seed = parseInt(this.getElement<HTMLInputElement>('#de-seed').value, 10);
        const limit = parseInt(this.getElement<HTMLInputElement>('#de-limit').value, 10);

        return {
            groupby,
            group1,
            group2,
            tests,
            max_cells_per_group: Number.isFinite(maxCells) ? maxCells : undefined,
            seed: Number.isFinite(seed) ? seed : undefined,
            limit: Number.isFinite(limit) ? limit : undefined,
        };
    }

    private setError(message: string | null): void {
        const errorEl = this.getElement('#de-form-error');
        if (!message) {
            errorEl.innerHTML = '';
            return;
        }
        errorEl.innerHTML = `<div class="error-box">${this.escapeHtml(message)}</div>`;
    }

    private escapeHtml(text: string): string {
        return text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
    }

    private uniqSorted(items: string[]): string[] {
        return Array.from(new Set(items)).sort((a, b) => a.localeCompare(b));
    }

    /**
     * Reset the form to initial state
     */
    public reset(): void {
        this.allGroupValues = [];
        this.isSubmitting = false;

        this.getElement<HTMLInputElement>('#de-test-ttest').checked = true;
        this.getElement<HTMLInputElement>('#de-test-ranksum').checked = true;
        this.getElement<HTMLSelectElement>('#de-groupby').selectedIndex = 0;
        this.getElement<HTMLInputElement>('#de-value-filter').value = '';
        this.getElement<HTMLInputElement>('#de-one-vs-rest').checked = false;
        this.getElement<HTMLSelectElement>('#de-group1').innerHTML = '';
        this.getElement<HTMLSelectElement>('#de-group2').innerHTML = '';
        this.getElement<HTMLSelectElement>('#de-group1').disabled = true;
        this.getElement<HTMLSelectElement>('#de-group2').disabled = true;
        this.getElement<HTMLInputElement>('#de-max-cells').value = '2000';
        this.getElement<HTMLInputElement>('#de-seed').value = '0';
        this.getElement<HTMLInputElement>('#de-limit').value = '50';

        this.updateChips();
        this.updateSubmitState();
        this.setError(null);
    }
}
