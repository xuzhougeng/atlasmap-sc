// SOMA-Tiles Frontend Entry Point

import { MapController } from './map/MapController';
import { GeneSelector } from './components/GeneSelector';
import { CategoryFilter } from './components/CategoryFilter';
import { StateManager, AppState } from './state/StateManager';
import { ApiClient } from './api/client';

// Initialize application
async function init() {
    console.log('Initializing SOMA-Tiles...');

    // Create API client
    const api = new ApiClient('/api');

    // Load metadata
    let metadata;
    try {
        metadata = await api.getMetadata();
        console.log('Loaded metadata:', metadata);
        updateDatasetInfo(metadata);
    } catch (error) {
        console.error('Failed to load metadata:', error);
        showError('Failed to connect to server');
        return;
    }

    // Initialize state manager
    const initialState: AppState = {
        zoom: 2,
        center: [128, 128],
        selectedCategories: [],
        expressionFilter: { gene: null, min: 0, max: 1 },
        selection: { id: null, type: null, polygon: null, cellCount: 0 },
        colorMode: 'category',
        colorGene: null,
        colorScale: 'viridis',
    };
    const state = new StateManager(initialState);

    // Initialize map
    const mapContainer = document.getElementById('map-container')!;
    const mapController = new MapController(mapContainer, {
        apiUrl: '',
        tileSize: 256,
        maxZoom: metadata.zoom_levels - 1,
        bounds: metadata.bounds,
    });

    // Initialize gene selector
    const geneSelector = new GeneSelector(
        document.getElementById('gene-selector')!,
        api
    );
    geneSelector.onGeneSelect((gene) => {
        console.log('Gene selected:', gene);
        state.setState({ colorMode: 'expression', colorGene: gene });
        mapController.setExpressionGene(gene);
    });

    // Initialize category filter
    if (metadata.categories) {
        const categoryFilter = new CategoryFilter(
            document.getElementById('category-filters')!,
            metadata.categories
        );
        categoryFilter.onFilterChange((categories) => {
            console.log('Categories filtered:', categories);
            state.setState({ selectedCategories: categories });
            // Refresh tiles with filter
        });
    }

    // Set up toolbar buttons
    setupToolbar(mapController, state);

    console.log('SOMA-Tiles initialized successfully');
}

function updateDatasetInfo(metadata: any) {
    const infoEl = document.getElementById('dataset-info');
    if (infoEl) {
        infoEl.innerHTML = `
            <span class="info-item">${metadata.dataset_name || 'Dataset'}</span>
            <span class="info-item">${formatNumber(metadata.n_cells)} cells</span>
            <span class="info-item">${metadata.n_genes_preaggregated} genes</span>
        `;
    }
}

function setupToolbar(map: MapController, state: StateManager) {
    const lassoBtn = document.getElementById('btn-lasso');
    const resetBtn = document.getElementById('btn-reset');

    if (lassoBtn) {
        lassoBtn.addEventListener('click', () => {
            // Toggle lasso mode
            lassoBtn.classList.toggle('active');
            // map.toggleLassoMode();
        });
    }

    if (resetBtn) {
        resetBtn.addEventListener('click', () => {
            map.resetView();
            state.setState({
                colorMode: 'category',
                colorGene: null,
                selectedCategories: [],
                selection: { id: null, type: null, polygon: null, cellCount: 0 },
            });
        });
    }
}

function showError(message: string) {
    const app = document.getElementById('app');
    if (app) {
        const error = document.createElement('div');
        error.className = 'error-message';
        error.textContent = message;
        app.prepend(error);
    }
}

function formatNumber(n: number): string {
    if (n >= 1e6) return (n / 1e6).toFixed(1) + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(1) + 'K';
    return n.toString();
}

// Start application
document.addEventListener('DOMContentLoaded', init);
