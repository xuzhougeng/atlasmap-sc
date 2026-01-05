// Differential Expression Page Entry Point

import './styles/main.css';

import { ApiClient } from './api/client';
import { ThemeManager } from './components/ThemeManager';
import { DEJobForm } from './components/DEJobForm';

interface DatasetInfo {
    id: string;
    name: string;
}

interface DatasetsResponse {
    default: string;
    datasets: DatasetInfo[];
    title: string;
}

function getDatasetFromUrl(): string | null {
    const params = new URLSearchParams(window.location.search);
    return params.get('dataset');
}

async function fetchDatasets(): Promise<DatasetsResponse> {
    const response = await fetch('/api/datasets');
    if (!response.ok) {
        throw new Error(`Failed to fetch datasets: ${response.statusText}`);
    }
    return response.json();
}

function createDatasetSelector(
    container: HTMLElement,
    datasets: DatasetInfo[],
    currentDataset: string
): void {
    if (datasets.length <= 1) return;

    const selectorContainer = document.createElement('div');
    selectorContainer.className = 'dataset-selector';

    const label = document.createElement('label');
    label.textContent = 'Dataset: ';
    label.htmlFor = 'dataset-select';

    const select = document.createElement('select');
    select.id = 'dataset-select';
    select.className = 'dataset-select';

    for (const ds of datasets) {
        const option = document.createElement('option');
        option.value = ds.id;
        option.textContent = ds.name;
        option.selected = ds.id === currentDataset;
        select.appendChild(option);
    }

    select.addEventListener('change', () => {
        const url = new URL(window.location.href);
        url.searchParams.set('dataset', select.value);
        window.location.href = url.toString();
    });

    selectorContainer.appendChild(label);
    selectorContainer.appendChild(select);
    container.appendChild(selectorContainer);
}

function showError(message: string): void {
    const root = document.getElementById('de-root');
    if (root) {
        root.innerHTML = `<div class="error-box">${message.replace(/</g, '&lt;')}</div>`;
    }
}

async function init(): Promise<void> {
    new ThemeManager();

    const root = document.getElementById('de-root');
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
        document.title = `${datasetsInfo.title} - DE`;
    }

    // Determine current dataset
    const urlDataset = getDatasetFromUrl();
    const validDatasetIds = datasetsInfo.datasets.map(d => d.id);
    const currentDataset = urlDataset && validDatasetIds.includes(urlDataset)
        ? urlDataset
        : datasetsInfo.default;

    // Setup dataset selector
    const datasetInfoContainer = document.getElementById('dataset-info');
    if (datasetInfoContainer) {
        createDatasetSelector(datasetInfoContainer, datasetsInfo.datasets, currentDataset);
    }

    // Update navigation links
    const backBtn = document.getElementById('btn-back') as HTMLAnchorElement | null;
    if (backBtn) {
        backBtn.href = `/?dataset=${encodeURIComponent(currentDataset)}`;
    }

    // Create API client and form
    const api = new ApiClient(`/d/${encodeURIComponent(currentDataset)}/api`);

    new DEJobForm(root, api, {
        onSubmit: (jobId) => {
            const url = new URL(window.location.href);
            url.pathname = '/de_job.html';
            url.searchParams.set('dataset', currentDataset);
            url.searchParams.set('job_id', jobId);
            window.location.href = url.toString();
        },
        onError: (error) => {
            console.error('DE job submission failed:', error);
        },
    });
}

document.addEventListener('DOMContentLoaded', () => {
    void init();
});
