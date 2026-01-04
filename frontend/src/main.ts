// SOMA-Tiles Frontend Entry Point

import { MapController } from './map/MapController';
import { GeneSelector } from './components/GeneSelector';
import { CategoryFilter } from './components/CategoryFilter';
import { TabPanel, TabType } from './components/TabPanel';
import { CategoryColumnSelector } from './components/CategoryColumnSelector';
import { CategoryLegend } from './components/CategoryLegend';
import { CellQueryPanel } from './components/CellQueryPanel';
import { ColorScaleSelector } from './components/ColorScaleSelector';
import { SidebarResizer } from './components/SidebarResizer';
import { StateManager, AppState } from './state/StateManager';
import { ApiClient } from './api/client';
import { downloadPngAtCurrentZoom } from './map/exportPng';

// Initialize application
async function init() {
    console.log('Initializing SOMA-Tiles...');

    // Create API client
    const api = new ApiClient('/api');
    let cellQueryPanel: CellQueryPanel | null = null;

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
    const maxNativeZoom = Math.max(0, metadata.zoom_levels - 1);
    // Allow a few extra zoom levels to "magnify" the highest native tiles.
    const maxZoom = maxNativeZoom + 4;
    const mapController = new MapController(mapContainer, {
        apiUrl: '',
        tileSize: 256,
        maxZoom,
        maxNativeZoom,
        // Start at an overview zoom so the full dataset is visible immediately.
        initialZoom: 0,
        bounds: metadata.bounds,
    });

    // Get available category columns
    const availableCategories = api.getAvailableCategories(metadata);
    const defaultCategory = availableCategories.includes('cell_type')
        ? 'cell_type'
        : availableCategories[0] || '';
    let currentCategoryColumn = defaultCategory;
    let categoryFilter: CategoryFilter | null = null;

    // Initialize TabPanel
    const tabPanel = new TabPanel(
        document.getElementById('tab-panel-container')!,
        {
            onTabChange: (tab: TabType) => {
                console.log('Tab changed:', tab);
                state.setState({ colorMode: tab });

                if (tab === 'category') {
                    const column = categoryColumnSelector.getSelectedCategory();
                    currentCategoryColumn = column;
                    const filter = categoryFilter?.getFilterForColumn(column) ?? null;
                    mapController.setCategoryColumn(column, filter);
                    categoryLegend.loadLegend(column);
                    categoryLegend.show();
                } else {
                    categoryLegend.hide();
                }
            },
        }
    );

    // Get tab containers
    const categoryContainer = tabPanel.getCategoryContainer();
    const expressionContainer = tabPanel.getExpressionContainer();

    // === Category Tab Content ===

    // Create category column selector container
    const categoryColumnSelectorContainer = document.createElement('div');
    categoryColumnSelectorContainer.id = 'category-column-selector';
    categoryContainer.appendChild(categoryColumnSelectorContainer);

    // Initialize category column selector
    const categoryColumnSelector = new CategoryColumnSelector(
        categoryColumnSelectorContainer,
        {
            onCategoryChange: (column) => {
                console.log('Category column changed:', column);
                currentCategoryColumn = column;
                const filter = categoryFilter?.getFilterForColumn(column) ?? null;
                mapController.setCategoryColumn(column, filter);
                categoryLegend.loadLegend(column);
                cellQueryPanel?.setCategoryColumn(column);
            },
        }
    );
    categoryColumnSelector.setCategories(availableCategories);
    if (defaultCategory) {
        categoryColumnSelector.setSelectedCategory(defaultCategory);
        currentCategoryColumn = defaultCategory;
    }

    // Create category legend container
    const categoryLegendContainer = document.createElement('div');
    categoryLegendContainer.id = 'category-legend';
    categoryLegendContainer.className = 'category-legend-container';
    categoryContainer.appendChild(categoryLegendContainer);

    // Initialize category legend
    const categoryLegend = new CategoryLegend(categoryLegendContainer, api);

    // Create category filter section
    const filterTitle = document.createElement('h3');
    filterTitle.className = 'tab-section-title';
    filterTitle.textContent = 'Filter';
    categoryContainer.appendChild(filterTitle);

    const categoryFiltersContainer = document.createElement('div');
    categoryFiltersContainer.id = 'category-filters';
    categoryFiltersContainer.className = 'category-filters';
    categoryContainer.appendChild(categoryFiltersContainer);

    // Initialize category filter
    if (metadata.categories) {
        categoryFilter = new CategoryFilter(
            categoryFiltersContainer,
            metadata.categories
        );
        categoryFilter.onFilterChange((column, filter) => {
            console.log('Categories filtered:', column, filter);
            if (column === currentCategoryColumn) {
                mapController.updateCategoryFilter(filter);
            }
            state.setState({ selectedCategories: filter ?? [] });
        });
    }

    // === Expression Tab Content ===

    // Create gene selector section
    const geneSelectorContainer = document.createElement('div');
    geneSelectorContainer.className = 'gene-selector';
    geneSelectorContainer.id = 'gene-selector';
    geneSelectorContainer.innerHTML = `
        <input type="text"
               id="gene-input"
               class="gene-input"
               placeholder="Search genes..."
               autocomplete="off">
        <div id="gene-dropdown" class="gene-dropdown"></div>
    `;
    expressionContainer.appendChild(geneSelectorContainer);

    // Initialize gene selector
    const geneSelector = new GeneSelector(geneSelectorContainer, api);

    // Create colormap selector container
    const colormapSelectorContainer = document.createElement('div');
    colormapSelectorContainer.id = 'colormap-selector';
    expressionContainer.appendChild(colormapSelectorContainer);

    // Initialize colormap selector
    const colormapSelector = new ColorScaleSelector(colormapSelectorContainer, {
        onScaleChange: (scale) => {
            state.setState({ colorScale: scale });
            const { colorGene, colorMode } = state.getState();
            if (colorMode === 'expression' && colorGene) {
                mapController.setExpressionGene(colorGene, scale);
            }
        },
    });
    colormapSelector.setScales(['viridis', 'plasma', 'inferno', 'magma']);
    colormapSelector.setSelectedScale(state.getState().colorScale);

    // Create expression legend container
    const expressionLegendContainer = document.createElement('div');
    expressionLegendContainer.id = 'expression-legend';
    expressionLegendContainer.className = 'legend';
    expressionContainer.appendChild(expressionLegendContainer);

    // Create cell query panel container
    const cellQueryContainer = document.createElement('div');
    cellQueryContainer.id = 'cell-query-panel';
    expressionContainer.appendChild(cellQueryContainer);

    // Initialize cell query panel
    cellQueryPanel = new CellQueryPanel(
        cellQueryContainer,
        api,
        {
            onBinSelect: (bin) => {
                console.log('Bin selected:', bin);
                mapController.panTo(bin.bin_x, bin.bin_y);
            },
        }
    );
    if (defaultCategory) {
        cellQueryPanel.setCategoryColumn(defaultCategory);
    }

    // Gene selection handler
    geneSelector.onGeneSelect((gene) => {
        console.log('Gene selected:', gene);
        state.setState({ colorMode: 'expression', colorGene: gene });

        // Switch to expression tab (this will also trigger mode change)
        tabPanel.switchTab('expression');
        mapController.setExpressionGene(gene, state.getState().colorScale);
        categoryLegend.hide();

        // Update cell query panel with selected gene
        cellQueryPanel?.setGene(gene);
    });

    // Initialize with default category
    if (defaultCategory) {
        const filter = categoryFilter?.getFilterForColumn(defaultCategory) ?? null;
        mapController.setCategoryColumn(defaultCategory, filter);
        categoryLegend.loadLegend(defaultCategory);
    }

    // Set up toolbar buttons
    setupToolbar(mapController, state, tabPanel, categoryLegend, cellQueryPanel!, defaultCategory);

    // Initialize sidebar resizer
    const sidebarResizer = document.getElementById('sidebar-resizer');
    const sidebar = document.querySelector('.sidebar') as HTMLElement;
    if (sidebarResizer && sidebar) {
        new SidebarResizer(sidebarResizer, sidebar);
    }

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

function setupToolbar(
    map: MapController,
    state: StateManager,
    tabPanel: TabPanel,
    categoryLegend: CategoryLegend,
    cellQueryPanel: CellQueryPanel,
    defaultCategory: string
) {
    const lassoBtn = document.getElementById('btn-lasso');
    const resetBtn = document.getElementById('btn-reset');
    const downloadBtn = document.getElementById('btn-download') as HTMLButtonElement | null;
    const downloadStatus = document.getElementById('download-status');

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

            // Reset to category tab
            tabPanel.switchTab('category');
            if (defaultCategory) {
                map.setCategoryColumn(defaultCategory);
                categoryLegend.loadLegend(defaultCategory);
                categoryLegend.show();
            }

            // Clear cell query panel
            cellQueryPanel.clear();

            state.setState({
                colorMode: 'category',
                colorGene: null,
                selectedCategories: [],
                selection: { id: null, type: null, polygon: null, cellCount: 0 },
            });
        });
    }

    if (downloadBtn) {
        let downloading = false;
        downloadBtn.addEventListener('click', async () => {
            if (downloading) return;
            downloading = true;
            downloadBtn.disabled = true;

            const mapZoom = map.getZoom();
            const exportZoom = Math.min(mapZoom, map.getMaxNativeZoom());
            if (downloadStatus) {
                downloadStatus.textContent =
                    mapZoom === exportZoom
                        ? `Preparing PNG (zoom ${exportZoom})...`
                        : `Preparing PNG (zoom ${exportZoom}, clamped from ${mapZoom})...`;
            }

            try {
                const result = await downloadPngAtCurrentZoom(map, {
                    maxDim: 8192,
                    concurrency: 12,
                    onProgress: ({ done, total }) => {
                        if (!downloadStatus) return;
                        downloadStatus.textContent = `Downloading tiles: ${done}/${total}`;
                    },
                });

                if (downloadStatus) {
                    downloadStatus.textContent = `Saved ${result.filename} (${result.width}×${result.height})`;
                    window.setTimeout(() => {
                        if (downloadStatus.textContent?.startsWith('Saved ')) {
                            downloadStatus.textContent = '';
                        }
                    }, 4000);
                }
            } catch (error) {
                console.error('Failed to download PNG:', error);
                if (downloadStatus) {
                    downloadStatus.textContent = 'Download failed (see console)';
                }
            } finally {
                downloading = false;
                downloadBtn.disabled = false;
            }
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
