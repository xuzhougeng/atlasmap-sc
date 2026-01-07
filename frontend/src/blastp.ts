// BLASTP Search Page Entry Point

import './styles/main.css';

import { ThemeManager } from './components/ThemeManager';

interface DatasetInfo {
    id: string;
    name: string;
    has_blastp?: boolean;
}

interface DatasetsResponse {
    default: string;
    datasets: DatasetInfo[];
    title: string;
}

async function fetchDatasets(): Promise<DatasetsResponse> {
    const response = await fetch('/api/datasets');
    if (!response.ok) {
        throw new Error(`Failed to fetch datasets: ${response.statusText}`);
    }
    return response.json();
}

function showError(message: string): void {
    const root = document.getElementById('blastp-root');
    if (root) {
        root.innerHTML = `<div class="error-box">${message.replace(/</g, '&lt;')}</div>`;
    }
}

function renderForm(root: HTMLElement, datasets: DatasetInfo[]): void {
    const blastpDatasets = datasets.filter(d => d.has_blastp);
    
    root.innerHTML = `
        <form id="blastp-form" class="blastp-form">
            <div class="blastp-form-group">
                <label class="blastp-label" for="sequence">Protein Sequence (FASTA or raw)</label>
                <textarea 
                    id="sequence" 
                    name="sequence" 
                    class="blastp-textarea" 
                    rows="8" 
                    placeholder=">query_protein
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQ..."
                    required
                ></textarea>
            </div>
            
            <div class="blastp-form-row">
                <div class="blastp-form-group">
                    <label class="blastp-label" for="max_hits">Max Hits per Database</label>
                    <input type="number" id="max_hits" name="max_hits" class="blastp-input" value="10" min="1" max="100" />
                </div>
                
                <div class="blastp-form-group">
                    <label class="blastp-label" for="evalue">E-value Cutoff</label>
                    <input type="text" id="evalue" name="evalue" class="blastp-input" value="1e-5" />
                </div>
            </div>
            
            ${blastpDatasets.length > 0 ? `
            <div class="blastp-form-group">
                <label class="blastp-label">Target Datasets (leave empty for all)</label>
                <div class="blastp-dataset-grid">
                    ${blastpDatasets.map(d => `
                        <label class="blastp-dataset-item">
                            <input type="checkbox" name="datasets" value="${d.id}" />
                            <span>${d.name}</span>
                        </label>
                    `).join('')}
                </div>
            </div>
            ` : `
            <div class="blastp-notice">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <circle cx="12" cy="12" r="10"/>
                    <line x1="12" y1="8" x2="12" y2="12"/>
                    <line x1="12" y1="16" x2="12.01" y2="16"/>
                </svg>
                <span>No datasets have BLASTP databases configured. Please configure blastp_path in server.yaml.</span>
            </div>
            `}
            
            <div class="blastp-form-actions">
                <button type="submit" class="blastp-submit-btn" ${blastpDatasets.length === 0 ? 'disabled' : ''}>
                    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <circle cx="11" cy="11" r="8"/>
                        <line x1="21" y1="21" x2="16.65" y2="16.65"/>
                    </svg>
                    Submit BLASTP Job
                </button>
            </div>
        </form>
    `;

    const form = document.getElementById('blastp-form') as HTMLFormElement;
    form?.addEventListener('submit', async (e) => {
        e.preventDefault();
        await submitJob(form);
    });
}

async function submitJob(form: HTMLFormElement): Promise<void> {
    const formData = new FormData(form);
    const sequence = formData.get('sequence') as string;
    const maxHits = parseInt(formData.get('max_hits') as string) || 10;
    const evalueStr = formData.get('evalue') as string;
    const evalue = parseFloat(evalueStr) || 1e-5;
    
    // Get selected datasets
    const selectedDatasets: string[] = [];
    form.querySelectorAll('input[name="datasets"]:checked').forEach((el) => {
        selectedDatasets.push((el as HTMLInputElement).value);
    });

    const submitBtn = form.querySelector('button[type="submit"]') as HTMLButtonElement;
    const originalText = submitBtn.innerHTML;
    submitBtn.disabled = true;
    submitBtn.innerHTML = 'Submitting...';

    try {
        const response = await fetch('/api/blastp/jobs', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                sequence,
                max_hits: maxHits,
                evalue,
                datasets: selectedDatasets.length > 0 ? selectedDatasets : undefined,
            }),
        });

        if (!response.ok) {
            const errText = await response.text();
            throw new Error(errText || response.statusText);
        }

        const result = await response.json();
        
        // Redirect to job status page
        const url = new URL(window.location.href);
        url.pathname = '/blastp_job.html';
        url.searchParams.set('job_id', result.job_id);
        window.location.href = url.toString();
    } catch (error: any) {
        console.error('Failed to submit BLASTP job:', error);
        alert('Failed to submit job: ' + (error?.message || error));
        submitBtn.disabled = false;
        submitBtn.innerHTML = originalText;
    }
}

async function init(): Promise<void> {
    new ThemeManager();

    const root = document.getElementById('blastp-root');
    if (!root) return;

    // Load datasets
    let datasetsInfo: DatasetsResponse;
    try {
        datasetsInfo = await fetchDatasets();
    } catch (e) {
        console.error(e);
        showError('Failed to connect to server (/api/datasets)');
        return;
    }

    // Update title
    if (datasetsInfo.title) {
        const titleEl = document.querySelector('.header-title');
        if (titleEl) titleEl.textContent = datasetsInfo.title;
        document.title = `${datasetsInfo.title} - BLASTP`;
    }

    // Render form
    renderForm(root, datasetsInfo.datasets);
}

document.addEventListener('DOMContentLoaded', () => {
    void init();
});

