// AtlasMap Frontend Entry Point

// Import CSS (Leaflet must come first so app styles can override)
import 'leaflet/dist/leaflet.css';
import './styles/main.css';

import { MapController } from './map/MapController';
import { CellCanvasLayer } from './map/CellCanvasLayer';
import { GeneSelector } from './components/GeneSelector';
import { CategoryFilter } from './components/CategoryFilter';
import { TabPanel, TabType } from './components/TabPanel';
import { CategoryColumnSelector } from './components/CategoryColumnSelector';
import { CategoryLegend } from './components/CategoryLegend';
import { CellQueryPanel } from './components/CellQueryPanel';
import { ColorScaleSelector } from './components/ColorScaleSelector';
import { ExpressionColorbar } from './components/ExpressionColorbar';
import { ExpressionRangeSelector } from './components/ExpressionRangeSelector';
import { SidebarResizer } from './components/SidebarResizer';
import { ThemeManager } from './components/ThemeManager';
import { LassoSelectionPanel } from './components/LassoSelectionPanel';
import { StateManager, AppState } from './state/StateManager';
import { ApiClient } from './api/client';
import { downloadPngAtCurrentZoom } from './map/exportPng';

interface DatasetInfo {
    id: string;
    name: string;
    has_h5ad?: boolean;
}

interface DatasetsResponse {
    default: string;
    datasets: DatasetInfo[];
    title: string;
}

// Flag to suppress labels toggle click when opening font-size popover via long-press
let labelLongPressJustHandled = false;

// Get current dataset from URL query params
function getCurrentDatasetFromUrl(): string | null {
    const params = new URLSearchParams(window.location.search);
    return params.get('dataset');
}

function getCurrentCoordFromUrl(): string | null {
    const params = new URLSearchParams(window.location.search);
    return params.get('coord');
}

// Set dataset in URL and reload
function setDatasetAndReload(datasetId: string): void {
    const url = new URL(window.location.href);
    url.searchParams.set('dataset', datasetId);
    window.location.href = url.toString();
}

function setCoordAndReload(coordKey: string): void {
    const url = new URL(window.location.href);
    url.searchParams.set('coord', coordKey);
    window.location.href = url.toString();
}

// Fetch available datasets from server
async function fetchDatasets(): Promise<DatasetsResponse> {
    const response = await fetch('/api/datasets');
    if (!response.ok) {
        throw new Error(`Failed to fetch datasets: ${response.statusText}`);
    }
    return response.json();
}

// Resolve gene_id to matching datasets
interface GeneLookupResponse {
    gene_id: string;
    datasets: string[];
}

async function lookupGene(geneId: string): Promise<GeneLookupResponse> {
    const response = await fetch(`/api/gene_lookup?gene_id=${encodeURIComponent(geneId)}`);
    if (!response.ok) {
        throw new Error(`Failed to lookup gene: ${response.statusText}`);
    }
    return response.json();
}

// Show gene dataset selection dialog when gene matches multiple datasets
function showGeneDatasetSelector(geneId: string, datasets: string[]): void {
    const overlay = document.createElement('div');
    overlay.className = 'gene-dataset-overlay';
    overlay.innerHTML = `
        <div class="gene-dataset-modal">
            <h2>Select Dataset</h2>
            <p>Gene <strong>${geneId}</strong> found in multiple datasets:</p>
            <div class="gene-dataset-list">
                ${datasets.map(ds => `
                    <button class="gene-dataset-btn" data-dataset="${ds}">${ds}</button>
                `).join('')}
            </div>
        </div>
    `;
    
    overlay.querySelectorAll('.gene-dataset-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            const dataset = (btn as HTMLElement).dataset.dataset;
            if (dataset) {
                const url = new URL(window.location.href);
                url.searchParams.set('dataset', dataset);
                url.searchParams.set('gene', geneId);
                window.location.href = url.toString();
            }
        });
    });
    
    document.body.appendChild(overlay);
}

// Create dataset selector dropdown
function createDatasetSelector(
    container: HTMLElement,
    datasets: DatasetInfo[],
    currentDataset: string
): HTMLDivElement | null {
    if (datasets.length <= 1) {
        // Only one dataset, no need for selector
        return null;
    }

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
        if (ds.id === currentDataset) {
            option.selected = true;
        }
        select.appendChild(option);
    }

    select.addEventListener('change', () => {
        setDatasetAndReload(select.value);
    });

    selectorContainer.appendChild(label);
    selectorContainer.appendChild(select);

    // Insert before the existing dataset-info content
    container.insertBefore(selectorContainer, container.firstChild);
    return selectorContainer;
}

function createCoordSelector(
    container: HTMLElement,
    coords: { key: string }[] | undefined,
    currentCoord: string
): HTMLDivElement | null {
    if (!coords || coords.length <= 1) {
        return null;
    }

    const selectorContainer = document.createElement('div');
    selectorContainer.className = 'dataset-selector';

    const label = document.createElement('label');
    label.textContent = 'Coord: ';
    label.htmlFor = 'coord-select';

    const select = document.createElement('select');
    select.id = 'coord-select';
    select.className = 'dataset-select';

    for (const cs of coords) {
        const option = document.createElement('option');
        option.value = cs.key;
        option.textContent = cs.key;
        option.selected = cs.key === currentCoord;
        select.appendChild(option);
    }

    select.addEventListener('change', () => {
        setCoordAndReload(select.value);
    });

    selectorContainer.appendChild(label);
    selectorContainer.appendChild(select);

    const selectors = container.querySelectorAll('.dataset-selector');
    const anchor = selectors.length ? selectors[selectors.length - 1] : null;
    if (anchor && anchor.parentElement === container) {
        container.insertBefore(selectorContainer, anchor.nextSibling);
    } else {
        container.insertBefore(selectorContainer, container.firstChild);
    }

    return selectorContainer;
}

function createH5ADDownloadButton(
    container: HTMLElement,
    datasetId: string,
    hasH5AD: boolean
): void {
    const existing = container.querySelector('#btn-download-h5ad');
    if (existing) existing.remove();

    if (!hasH5AD) return;

    const link = document.createElement('a');
    link.id = 'btn-download-h5ad';
    link.className = 'tool-btn header-h5ad-download-btn';
    link.href = `/d/${encodeURIComponent(datasetId)}/api/h5ad`;
    link.download = `${datasetId}.h5ad`;
    link.title = 'Download H5AD';
    link.setAttribute('aria-label', 'Download H5AD');
    link.innerHTML = `
        <svg viewBox="0 0 24 24" width="20" height="20" aria-hidden="true" focusable="false">
            <path fill="currentColor" d="M5 20h14v-2H5v2zM12 2v12l4-4 1.41 1.41L12 17.83 6.59 11.41 8 10l4 4V2h0z"/>
        </svg>
    `;

    const selectors = container.querySelectorAll('.dataset-selector');
    const anchor = selectors.length ? selectors[selectors.length - 1] : null;
    if (anchor && anchor.parentElement === container) {
        container.insertBefore(link, anchor.nextSibling);
    } else {
        container.appendChild(link);
    }
}

// Initialize theme manager early to prevent flash
const themeManager = new ThemeManager();

function setupThemeToggle() {
    const themeBtn = document.getElementById('btn-theme');
    const iconSun = document.getElementById('icon-sun');
    const iconMoon = document.getElementById('icon-moon');

    const updateIcon = (theme: 'dark' | 'light') => {
        if (iconSun && iconMoon) {
            iconSun.style.display = theme === 'dark' ? 'block' : 'none';
            iconMoon.style.display = theme === 'light' ? 'block' : 'none';
        }
    };

    updateIcon(themeManager.getTheme());

    themeBtn?.addEventListener('click', () => {
        const newTheme = themeManager.toggle();
        updateIcon(newTheme);
    });
}

// Initialize application
async function init() {
    console.log('Initializing AtlasMap...');
    setupThemeToggle();

    // Fetch available datasets
    let datasetsInfo: DatasetsResponse;
    try {
        datasetsInfo = await fetchDatasets();
        console.log('Available datasets:', datasetsInfo);
    } catch (error) {
        console.error('Failed to fetch datasets:', error);
        showError('Failed to connect to server');
        return;
    }

    // Update page title from config
    if (datasetsInfo.title) {
        const titleEl = document.querySelector('.header-title');
        if (titleEl) {
            titleEl.textContent = datasetsInfo.title;
        }
        document.title = datasetsInfo.title;
    }

    // Check for gene auto-resolution: if gene is provided but dataset is missing/invalid
    const urlParams = new URLSearchParams(window.location.search);
    const geneFromUrl = urlParams.get('gene');
    const urlDataset = getCurrentDatasetFromUrl();
    const validDatasetIds = datasetsInfo.datasets.map(d => d.id);

    // If gene is provided but dataset is missing or invalid, try to auto-resolve
    if (geneFromUrl && (!urlDataset || !validDatasetIds.includes(urlDataset))) {
        try {
            console.log('Auto-resolving dataset for gene:', geneFromUrl);
            const lookup = await lookupGene(geneFromUrl);
            
            if (lookup.datasets.length === 1) {
                // Single match: redirect automatically
                const url = new URL(window.location.href);
                url.searchParams.set('dataset', lookup.datasets[0]);
                url.searchParams.set('gene', geneFromUrl);
                window.location.href = url.toString();
                return;
            } else if (lookup.datasets.length > 1) {
                // Multiple matches: show selection dialog
                showGeneDatasetSelector(geneFromUrl, lookup.datasets);
                return;
            } else {
                // No match: show warning but continue with default dataset
                console.warn(`Gene ${geneFromUrl} not found in any dataset`);
            }
        } catch (error) {
            console.error('Failed to lookup gene:', error);
            // Continue with normal flow
        }
    }

    // Determine which dataset to use
    let currentDataset: string;

    if (urlDataset && validDatasetIds.includes(urlDataset)) {
        currentDataset = urlDataset;
    } else {
        currentDataset = datasetsInfo.default;
        // If URL had an invalid dataset, redirect to the default
        if (urlDataset && urlDataset !== currentDataset) {
            setDatasetAndReload(currentDataset);
            return;
        }
    }

    console.log('Using dataset:', currentDataset);

    // Compute API base URL for the selected dataset
    const apiBaseUrl = `/d/${encodeURIComponent(currentDataset)}/api`;
    const tilesBaseUrl = `/d/${encodeURIComponent(currentDataset)}`;

    // Create API client
    const urlCoord = getCurrentCoordFromUrl();
    let api = new ApiClient(apiBaseUrl, urlCoord ? { coord: urlCoord } : {});
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

    const effectiveCoord =
        metadata.coordinate_system ||
        urlCoord ||
        metadata.default_coordinate_system ||
        metadata.umap_key ||
        null;

    // Keep the URL in sync if the server fell back to a different coord.
    if (effectiveCoord && urlCoord && effectiveCoord !== urlCoord) {
        const url = new URL(window.location.href);
        url.searchParams.set('coord', effectiveCoord);
        window.history.replaceState({}, '', url.toString());
    }

    // Recreate API client so all subsequent API calls use the effective coord.
    if (effectiveCoord) {
        api = new ApiClient(apiBaseUrl, { coord: effectiveCoord });
    }

    // Create dataset selector in header
    const datasetInfoContainer = document.getElementById('dataset-info');
    if (datasetInfoContainer) {
        createDatasetSelector(datasetInfoContainer, datasetsInfo.datasets, currentDataset);
        createCoordSelector(datasetInfoContainer, metadata.coordinate_systems, effectiveCoord ?? '');
        const currentDatasetInfo = datasetsInfo.datasets.find(ds => ds.id === currentDataset);
        createH5ADDownloadButton(datasetInfoContainer, currentDataset, currentDatasetInfo?.has_h5ad === true);
    }

    // Initialize state manager
    const initialState: AppState = {
        zoom: 2,
        center: [128, 128],
        selectedCategories: [],
        expressionFilter: { gene: null, min: 0, max: 1 },
        colorMode: 'category',
        colorGene: null,
        colorScale: 'seurat',
    };
    const state = new StateManager(initialState);

    // Initialize map
    const mapContainer = document.getElementById('map-container')!;
    
    // renderZoom is the maximum bins resolution from preprocessing (zoom_levels - 1)
    const renderZoom = Math.max(0, metadata.zoom_levels - 1);
    
    const tileSize = 256;
    // Base zoom limit for tiles (0..4 by default, clamped to renderZoom).
    const baseMaxNativeZoom = Math.min(4, renderZoom);
    const canvasExtraZoom = 4;

    let currentMaxNativeZoom = baseMaxNativeZoom;
    let currentMaxZoom = baseMaxNativeZoom;

    const mapController = new MapController(mapContainer, {
        apiUrl: tilesBaseUrl,
        coord: effectiveCoord ?? undefined,
        tileSize,
        renderZoom,
        maxZoom: currentMaxZoom,
        maxNativeZoom: currentMaxNativeZoom,
        // Start at zoom level 2 for better initial detail.
        initialZoom: 2,
        bounds: metadata.bounds,
    });

    // Initialize cell canvas layer for rendering at zoom > maxNativeZoom
    const cellCanvasLayer = new CellCanvasLayer({
        api,
        map: mapController.getMap(),
        maxNativeZoom: currentMaxNativeZoom,
        bounds: metadata.bounds,
    });
    // Default: cell canvas disabled
    cellCanvasLayer.setEnabled(false);
    
    // Apply cell canvas toggle: updates map constraints and cell canvas layer.
    function applyCellCanvasMode(enabled: boolean): void {
        const newMaxNativeZoom = baseMaxNativeZoom;
        const newMaxZoom = enabled ? baseMaxNativeZoom + canvasExtraZoom : baseMaxNativeZoom;

        currentMaxNativeZoom = newMaxNativeZoom;
        currentMaxZoom = newMaxZoom;

        // Update map zoom constraints (this also rebuilds tile layers)
        mapController.setZoomConstraints(newMaxZoom, newMaxNativeZoom);

        // Update cell canvas layer
        cellCanvasLayer.setEnabled(enabled);
        cellCanvasLayer.setMaxNativeZoom(newMaxNativeZoom);

        console.log(`Cell canvas: ${enabled ? 'ON' : 'OFF'}, maxNativeZoom=${newMaxNativeZoom}, maxZoom=${newMaxZoom}`);
    }

    // Helper to update cell canvas layer for category mode
    async function updateCellCanvasLayerCategory(column: string, filter: string[] | null): Promise<void> {
        try {
            const colors = await api.getCategoryColors(column);
            cellCanvasLayer.setColorOptions({
                mode: 'category',
                categoryColumn: column,
                categoryFilter: filter,
                categoryColors: colors,
            });
        } catch (error) {
            console.warn('Failed to load category colors for cell layer:', error);
            cellCanvasLayer.setColorOptions({
                mode: 'category',
                categoryColumn: column,
                categoryFilter: filter,
            });
        }
    }

    // Helper to update cell canvas layer for expression mode
    async function updateCellCanvasLayerExpression(gene: string, colormap: string): Promise<void> {
        try {
            // Use mapController.getMaxNativeZoom() to get the current value (may change with zoom mode)
            const stats = await api.getGeneStats(gene, mapController.getMaxNativeZoom());
            const p80 = stats.p80_expression;
            const autoMax = typeof p80 === 'number' && Number.isFinite(p80) && p80 > 0 ? p80 : stats.max_expression;
            cellCanvasLayer.setColorOptions({
                mode: 'expression',
                gene,
                expressionColormap: colormap,
                expressionMin: 0,
                expressionMax: autoMax,
            });
        } catch (error) {
            console.warn('Failed to load gene stats for cell layer:', error);
            cellCanvasLayer.setColorOptions({
                mode: 'expression',
                gene,
                expressionColormap: colormap,
            });
        }
    }

    // Expression colorbar overlay (shown only when expression tiles are active)
    const expressionColorbarContainer = document.createElement('div');
    expressionColorbarContainer.id = 'expression-colorbar';
    mapContainer.appendChild(expressionColorbarContainer);
    const expressionColorbar = new ExpressionColorbar(expressionColorbarContainer);

    let expressionColorbarStatsRequestId = 0;
    let expressionColorbarAutoCache: { gene: string; zoom: number; max: number } | null = null;

    // Get available category columns
    const availableCategories = api.getAvailableCategories(metadata);
    const defaultCategory = availableCategories.includes('cell_type')
        ? 'cell_type'
        : availableCategories[0] || '';
    let currentCategoryColumn = defaultCategory;
    let categoryFilter: CategoryFilter | null = null;
    let currentCategoryFilter: string[] | null = null;
    let showCategoryLabels = true;

    const labelsBtn = document.getElementById('btn-labels') as HTMLButtonElement | null;
    if (labelsBtn) {
        labelsBtn.classList.toggle('active', showCategoryLabels);
        labelsBtn.addEventListener('click', () => {
            if (labelLongPressJustHandled) {
                labelLongPressJustHandled = false;
                return;
            }
            showCategoryLabels = !showCategoryLabels;
            labelsBtn.classList.toggle('active', showCategoryLabels);
            void refreshCategoryLabels();
        });
    }

    // Cell canvas toggle button (OFF by default)
    const zoomSwitchBtn = document.getElementById('btn-zoom-switch') as HTMLButtonElement | null;
    const iconZoomTiles = document.getElementById('icon-zoom-tiles') as HTMLElement | null;
    const iconZoomCells = document.getElementById('icon-zoom-cells') as HTMLElement | null;
    let cellCanvasEnabled = false;

    function updateCellCanvasButton(enabled: boolean): void {
        if (iconZoomTiles && iconZoomCells) {
            iconZoomTiles.style.display = enabled ? 'none' : 'block';
            iconZoomCells.style.display = enabled ? 'block' : 'none';
        }
        if (zoomSwitchBtn) {
            zoomSwitchBtn.title = enabled
                ? 'Cell canvas: ON (allow zoom beyond tiles; tiles scale)'
                : `Cell canvas: OFF (zoom locked to tiles 0-${baseMaxNativeZoom})`;
            zoomSwitchBtn.classList.toggle('active', enabled);
        }
    }

    if (zoomSwitchBtn) {
        updateCellCanvasButton(cellCanvasEnabled);
        zoomSwitchBtn.addEventListener('click', () => {
            cellCanvasEnabled = !cellCanvasEnabled;
            applyCellCanvasMode(cellCanvasEnabled);
            updateCellCanvasButton(cellCanvasEnabled);

            // Update cellQueryPanel if it exists (defined later in init)
            if (typeof cellQueryPanel !== 'undefined' && cellQueryPanel) {
                cellQueryPanel.setMaxNativeZoom(mapController.getMaxNativeZoom());
            }
        });
    }

    // Prompt to enable cell canvas when user tries to zoom past tile limit.
    const zoomPrompt = document.createElement('div');
    zoomPrompt.className = 'zoom-canvas-prompt hidden';
    const zoomRangeLabel = `0-${baseMaxNativeZoom}`;
    zoomPrompt.innerHTML =
        `<div class="zoom-canvas-prompt__text">Reached max tile zoom (${zoomRangeLabel}). Enable cell canvas to zoom further?</div>` +
        '<div class="zoom-canvas-prompt__actions">' +
        '<button type="button" class="zoom-canvas-prompt__btn zoom-canvas-prompt__btn--primary">Enable</button>' +
        '<button type="button" class="zoom-canvas-prompt__btn zoom-canvas-prompt__btn--ghost">Dismiss</button>' +
        '</div>';
    mapContainer.appendChild(zoomPrompt);

    const enableBtn = zoomPrompt.querySelector('.zoom-canvas-prompt__btn--primary') as HTMLButtonElement | null;
    const dismissBtn = zoomPrompt.querySelector('.zoom-canvas-prompt__btn--ghost') as HTMLButtonElement | null;
    let zoomPromptTimer: number | null = null;

    const hideZoomPrompt = (): void => {
        zoomPrompt.classList.add('hidden');
        if (zoomPromptTimer !== null) {
            clearTimeout(zoomPromptTimer);
            zoomPromptTimer = null;
        }
    };

    const showZoomPrompt = (): void => {
        if (cellCanvasEnabled) return;
        if (window.matchMedia('(max-width: 768px)').matches) return; // skip on mobile
        zoomPrompt.classList.remove('hidden');
        if (zoomPromptTimer !== null) {
            clearTimeout(zoomPromptTimer);
        }
        zoomPromptTimer = window.setTimeout(() => {
            zoomPromptTimer = null;
            zoomPrompt.classList.add('hidden');
        }, 3500);
    };

    const leafletMap = mapController.getMap();
    const tryPromptOnZoomIn = (): void => {
        if (cellCanvasEnabled) return;
        const currentZoom = leafletMap.getZoom();
        const maxZoom = leafletMap.getMaxZoom();
        if (currentZoom >= maxZoom) {
            showZoomPrompt();
        }
    };

    mapContainer.addEventListener('wheel', (event) => {
        if (event.deltaY < 0) {
            tryPromptOnZoomIn();
        }
    }, { passive: true });

    mapContainer.addEventListener('dblclick', () => {
        tryPromptOnZoomIn();
    });

    const zoomInControl = mapContainer.querySelector('.leaflet-control-zoom-in') as HTMLAnchorElement | null;
    if (zoomInControl) {
        zoomInControl.addEventListener('click', () => {
            tryPromptOnZoomIn();
        });
    }

    enableBtn?.addEventListener('click', () => {
        if (!cellCanvasEnabled) {
            cellCanvasEnabled = true;
            applyCellCanvasMode(true);
            updateCellCanvasButton(true);
            if (typeof cellQueryPanel !== 'undefined' && cellQueryPanel) {
                cellQueryPanel.setMaxNativeZoom(mapController.getMaxNativeZoom());
            }
        }
        hideZoomPrompt();
        leafletMap.zoomIn(1);
    });

    dismissBtn?.addEventListener('click', () => {
        hideZoomPrompt();
    });

    leafletMap.on('zoomend', () => {
        hideZoomPrompt();
    });

    const refreshCategoryLabels = async (): Promise<void> => {
        mapController.setCategoryLabelsVisible(showCategoryLabels);
        if (!showCategoryLabels) {
            mapController.setCategoryLabelItems(null);
            return;
        }
        if (!currentCategoryColumn) {
            mapController.setCategoryLabelItems(null);
            return;
        }
        try {
            const items = await api.getCategoryCentroids(currentCategoryColumn);
            mapController.setCategoryLabelItems(items);
            mapController.setCategoryLabelFilter(currentCategoryFilter);
        } catch (error) {
            console.warn('Failed to load category centroids:', error);
            mapController.setCategoryLabelItems(null);
        }
    };

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
                    categoryFilter?.setActiveColumn(column);
                    const filter = categoryFilter?.getFilterForColumn(column) ?? null;
                    mapController.setCategoryColumn(column, filter);
                    categoryLegend.loadLegend(column);
                    categoryLegend.show();
                    // Update cell canvas layer for category mode
                    void updateCellCanvasLayerCategory(column, filter);
                } else {
                    categoryLegend.hide();
                    // Cell canvas will be updated when gene is selected
                }

                void refreshCategoryLabels();
                void refreshExpressionColorbar();
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
                categoryFilter?.setActiveColumn(column);
                const filter = categoryFilter?.getFilterForColumn(column) ?? null;
                mapController.setCategoryColumn(column, filter);
                categoryLegend.loadLegend(column);
                cellQueryPanel?.setCategoryColumn(column);
                currentCategoryFilter = filter ?? null;
                void refreshCategoryLabels();
                // Update cell canvas layer for category mode
                void updateCellCanvasLayerCategory(column, filter);
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
            metadata.categories,
            currentCategoryColumn
        );
        categoryFilter.onFilterChange((column, filter) => {
            console.log('Categories filtered:', column, filter);
            if (column === currentCategoryColumn) {
                mapController.updateCategoryFilter(filter);
                currentCategoryFilter = filter ?? null;
                mapController.setCategoryLabelFilter(currentCategoryFilter);
                // Update cell canvas layer filter
                cellCanvasLayer.setColorOptions({ categoryFilter: filter ?? null });
            }
            if (column === currentCategoryColumn) {
                state.setState({ selectedCategories: filter ?? [] });
            }
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

    let expressionRangeSelector: ExpressionRangeSelector | null = null;

    async function refreshExpressionColorbar(): Promise<void> {
        const mode = mapController.getColorMode();
        const gene = mapController.getCurrentExpressionGene();
        if (mode !== 'expression' || !gene) {
            expressionColorbar.hide();
            return;
        }

        expressionColorbar.show();
        expressionColorbar.setGene(gene);
        expressionColorbar.setColormap(mapController.getCurrentExpressionColorScale() ?? state.getState().colorScale);

        const manualRange = expressionRangeSelector?.getRange() ?? null;
        if (manualRange) {
            expressionColorbar.setRange(manualRange.min, manualRange.max, 'manual');
            return;
        }

        // Auto range: show 0..p80 at the server's native max zoom (matches backend tile scaling).
        const statsZoom = mapController.getMaxNativeZoom();
        if (
            expressionColorbarAutoCache &&
            expressionColorbarAutoCache.gene === gene &&
            expressionColorbarAutoCache.zoom === statsZoom
        ) {
            expressionColorbar.setRange(0, expressionColorbarAutoCache.max, 'auto');
            return;
        }

            expressionColorbar.setRange(0, null, 'auto');
            const requestId = ++expressionColorbarStatsRequestId;
            try {
                const stats = await api.getGeneStats(gene, statsZoom);
                if (requestId !== expressionColorbarStatsRequestId) return;
            const p80 = stats.p80_expression;
            const autoMax = typeof p80 === 'number' && Number.isFinite(p80) && p80 > 0 ? p80 : stats.max_expression;
            expressionColorbarAutoCache = { gene, zoom: statsZoom, max: autoMax };
            expressionColorbar.setRange(0, autoMax, 'auto');
        } catch (error) {
            console.error('Failed to load gene stats for expression colorbar:', error);
            if (requestId !== expressionColorbarStatsRequestId) return;
            expressionColorbar.setRange(0, null, 'auto');
        }
    }

    // Initialize colormap selector
    const colormapSelector = new ColorScaleSelector(colormapSelectorContainer, {
        onScaleChange: (scale) => {
            state.setState({ colorScale: scale });
            const { colorGene, colorMode } = state.getState();
            if (colorMode === 'expression' && colorGene) {
                mapController.setExpressionGene(
                    colorGene,
                    scale,
                    expressionRangeSelector?.getRange() ?? null
                );
                // Update cell canvas layer colormap
                cellCanvasLayer.setColorOptions({ expressionColormap: scale });
            }
            void refreshExpressionColorbar();
        },
    });
    colormapSelector.setScales(['viridis', 'plasma', 'inferno', 'magma', 'seurat']);
    colormapSelector.setSelectedScale(state.getState().colorScale);

    // Create expression range selector container
    const expressionRangeSelectorContainer = document.createElement('div');
    expressionRangeSelectorContainer.id = 'expression-range-selector';
    expressionContainer.appendChild(expressionRangeSelectorContainer);

    // Initialize expression range selector
    expressionRangeSelector = new ExpressionRangeSelector(expressionRangeSelectorContainer, {
        onRangeChange: (range) => {
            const { colorGene, colorMode, colorScale } = state.getState();
            if (colorMode === 'expression' && colorGene) {
                mapController.setExpressionGene(colorGene, colorScale, range);
                // Update cell canvas layer range
                if (range) {
                    cellCanvasLayer.setColorOptions({
                        expressionMin: range.min,
                        expressionMax: range.max,
                    });
                }
            }
            void refreshExpressionColorbar();
        },
    });

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
    // Set max native zoom so stats are fetched at correct zoom levels
    cellQueryPanel.setMaxNativeZoom(mapController.getMaxNativeZoom());
    if (defaultCategory) {
        cellQueryPanel.setCategoryColumn(defaultCategory);
    }

    // Listen for zoom changes to update stats panel
    mapController.getMap().on('zoomend', () => {
        const zoom = mapController.getZoom();
        cellQueryPanel?.setZoom(zoom);
        void refreshExpressionColorbar();
    });

    // Gene selection handler
    geneSelector.onGeneSelect((gene) => {
        console.log('Gene selected:', gene);
        state.setState({ colorMode: 'expression', colorGene: gene });

        // Switch to expression tab (this will also trigger mode change)
        tabPanel.switchTab('expression');
        mapController.setExpressionGene(
            gene,
            state.getState().colorScale,
            expressionRangeSelector?.getRange() ?? null
        );
        categoryLegend.hide();
        void refreshExpressionColorbar();

        // Update cell query panel with selected gene
        cellQueryPanel?.setGene(gene);

        // Update cell canvas layer for expression mode
        void updateCellCanvasLayerExpression(gene, state.getState().colorScale);
    });

    // Handle gene parameter from URL (e.g., from DE results page)
    const geneUrlParams = new URLSearchParams(window.location.search);
    const geneParamFromUrl = geneUrlParams.get('gene');
    if (geneParamFromUrl) {
        console.log('Gene from URL:', geneParamFromUrl);
        state.setState({ colorMode: 'expression', colorGene: geneParamFromUrl });
        tabPanel.switchTab('expression');
        mapController.setExpressionGene(
            geneParamFromUrl,
            state.getState().colorScale,
            expressionRangeSelector?.getRange() ?? null
        );
        categoryLegend.hide();
        void refreshExpressionColorbar();
        cellQueryPanel?.setGene(geneParamFromUrl);

        // Update gene input to show the selected gene
        const geneInput = document.getElementById('gene-input') as HTMLInputElement | null;
        if (geneInput) {
            geneInput.value = geneParamFromUrl;
        }

        // Update cell canvas layer for expression mode
        void updateCellCanvasLayerExpression(geneParamFromUrl, state.getState().colorScale);
    }

    // Initialize with default category (skip if gene was specified in URL)
    if (defaultCategory && !geneFromUrl) {
        const filter = categoryFilter?.getFilterForColumn(defaultCategory) ?? null;
        mapController.setCategoryColumn(defaultCategory, filter);
        categoryLegend.loadLegend(defaultCategory);
        currentCategoryFilter = filter ?? null;
        void refreshCategoryLabels();
        // Initialize cell canvas layer for category mode
        void updateCellCanvasLayerCategory(defaultCategory, filter);
    }

    // Set up toolbar buttons
    setupToolbar(
        mapController,
        state,
        tabPanel,
        categoryLegend,
        cellQueryPanel!,
        defaultCategory,
        currentDataset,
        () => {
            currentCategoryColumn = defaultCategory;
            currentCategoryFilter = null;
            mapController.setCategoryLabelFilter(null);
            void refreshCategoryLabels();
            void refreshExpressionColorbar();
        }
    );

    // Lasso selection tool + selection summary panel
    const selectionPanelContent = document.getElementById('selection-panel-content');
    const lassoBtn = document.getElementById('btn-lasso') as HTMLButtonElement | null;
    if (selectionPanelContent) {
        const lassoPanel = new LassoSelectionPanel(selectionPanelContent, mapController, api, {
            getCategoryColumn: () => currentCategoryColumn,
            getCategoryFilter: () => currentCategoryFilter,
            onEnabledChange: (enabled) => {
                if (!lassoBtn) return;
                lassoBtn.classList.toggle('active', enabled);
                lassoBtn.title = enabled ? 'Lasso select: ON (Esc to exit)' : 'Lasso select (lock pan/zoom)';
            },
        });

        if (lassoBtn) {
            lassoBtn.classList.toggle('active', lassoPanel.isEnabled());
            lassoBtn.addEventListener('click', () => lassoPanel.toggle());
        }
    } else if (lassoBtn) {
        lassoBtn.disabled = true;
        lassoBtn.title = 'Lasso select unavailable (missing selection panel)';
    }

    // Initialize sidebar resizer
    const sidebarResizer = document.getElementById('sidebar-resizer');
    const sidebar = document.querySelector('.sidebar') as HTMLElement;
    if (sidebarResizer && sidebar) {
        new SidebarResizer(sidebarResizer, sidebar);
    }

    // Handle viewport resize for mobile (debounced via requestAnimationFrame)
    let resizeRafId: number | null = null;
    const handleResize = () => {
        if (resizeRafId !== null) return;
        resizeRafId = requestAnimationFrame(() => {
            resizeRafId = null;
            mapController.getMap().invalidateSize({ pan: false });
        });
    };
    window.addEventListener('resize', handleResize);

    // Also handle orientation change on mobile devices
    window.addEventListener('orientationchange', () => {
        // Delay invalidateSize to let the browser settle after orientation change
        setTimeout(() => {
            mapController.getMap().invalidateSize({ pan: false });
        }, 100);
    });

    // Initial invalidateSize after layout paint (helps on mobile Safari)
    requestAnimationFrame(() => {
        requestAnimationFrame(() => {
            mapController.getMap().invalidateSize({ pan: false });
        });
    });

    // Setup mobile floating UI (only activates on mobile viewports)
    setupMobileUI(tabPanel, mapController);

    console.log('AtlasMap initialized successfully');
}

function updateDatasetInfo(metadata: any) {
    const infoEl = document.getElementById('dataset-info');
    if (infoEl) {
        const infoSpans = `
            <span class="info-item">${metadata.dataset_name || 'Dataset'}</span>
            <span class="info-item">${formatNumber(metadata.n_cells)} cells</span>
            <span class="info-item">${metadata.n_genes_preaggregated} genes</span>
        `;
        let infoContainer = infoEl.querySelector('.dataset-info-items') as HTMLSpanElement | null;
        if (!infoContainer) {
            infoContainer = document.createElement('span');
            infoContainer.className = 'dataset-info-items';
            infoEl.appendChild(infoContainer);
        }
        infoContainer.innerHTML = infoSpans;
    }
}

function setupToolbar(
    map: MapController,
    state: StateManager,
    tabPanel: TabPanel,
    categoryLegend: CategoryLegend,
    cellQueryPanel: CellQueryPanel,
    defaultCategory: string,
    datasetId: string,
    onAfterReset?: () => void
) {
    const resetBtn = document.getElementById('btn-reset');
    const panLockBtn = document.getElementById('btn-pan-lock') as HTMLButtonElement | null;
    const panLockedIcon = document.getElementById('icon-pan-locked') as HTMLElement | null;
    const panUnlockedIcon = document.getElementById('icon-pan-unlocked') as HTMLElement | null;
    const downloadBtn = document.getElementById('btn-download') as HTMLButtonElement | null;
    const deBtn = document.getElementById('btn-de') as HTMLButtonElement | null;
    const downloadStatus = document.getElementById('download-status');
    const pointSizeBtn = document.getElementById('btn-point-size') as HTMLButtonElement | null;
    const pointSizeSlider = document.getElementById('point-size-slider') as HTMLInputElement | null;
    const pointSizeValue = document.getElementById('point-size-value');
    const pointSizePopover = document.getElementById('point-size-popover') as HTMLDivElement | null;
    const labelsBtn = document.getElementById('btn-labels') as HTMLButtonElement | null;
    const labelFontSizePopover = document.getElementById('label-font-size-popover') as HTMLDivElement | null;
    const labelFontSizeSlider = document.getElementById('label-font-size-slider') as HTMLInputElement | null;
    const labelFontSizeValue = document.getElementById('label-font-size-value');

    const updatePanLockUi = () => {
        if (!panLockBtn) return;
        const locked = map.isPanLocked();
        panLockBtn.classList.toggle('active', !locked);
        panLockBtn.title = locked ? 'Unlock panning (allow free drag)' : 'Lock panning (stay within data bounds)';
        if (panLockedIcon && panUnlockedIcon) {
            panLockedIcon.style.display = locked ? 'block' : 'none';
            panUnlockedIcon.style.display = locked ? 'none' : 'block';
        }
    };

    const setPointSizeUi = (size: number) => {
        if (!pointSizeValue) return;
        pointSizeValue.textContent = `${Math.round(size * 100)}%`;
    };

    const applyPointSize = (size: number) => {
        map.setPointSize(size);
        const effective = map.getPointSize();
        if (pointSizeSlider) {
            pointSizeSlider.value = String(effective);
        }
        setPointSizeUi(effective);
    };

    if (pointSizeSlider) {
        const onInput = () => {
            const value = Number.parseFloat(pointSizeSlider.value);
            if (!Number.isFinite(value)) return;
            applyPointSize(value);
        };
        pointSizeSlider.addEventListener('input', onInput);
        onInput();
        pointSizeSlider.disabled = true;
    }

    const positionPointSizePopover = () => {
        if (!pointSizePopover || !pointSizeBtn) return;
        const panel = pointSizeBtn.closest('.panel') as HTMLElement | null;
        if (!panel) return;

        const panelRect = panel.getBoundingClientRect();
        const btnRect = pointSizeBtn.getBoundingClientRect();

        const padding = 8;
        const gap = 8;

        const popoverWidth = pointSizePopover.offsetWidth;
        const popoverHeight = pointSizePopover.offsetHeight;

        let left = btnRect.left - panelRect.left;
        let top = btnRect.bottom - panelRect.top + gap;

        const maxLeft = panelRect.width - popoverWidth - padding;
        if (left > maxLeft) left = maxLeft;
        if (left < padding) left = padding;

        if (top + popoverHeight > panelRect.height - padding) {
            top = btnRect.top - panelRect.top - popoverHeight - gap;
        }
        if (top < padding) top = padding;

        pointSizePopover.style.left = `${left}px`;
        pointSizePopover.style.top = `${top}px`;
    };

    const closePointSizePopover = () => {
        if (!pointSizePopover || !pointSizeBtn) return;
        pointSizePopover.classList.remove('open');
        pointSizePopover.setAttribute('aria-hidden', 'true');
        pointSizeBtn.classList.remove('active');
        if (pointSizeSlider) {
            pointSizeSlider.disabled = true;
        }
    };

    const openPointSizePopover = () => {
        if (!pointSizePopover || !pointSizeBtn) return;
        pointSizePopover.classList.add('open');
        pointSizePopover.setAttribute('aria-hidden', 'false');
        pointSizeBtn.classList.add('active');
        if (pointSizeSlider) {
            pointSizeSlider.disabled = false;
        }
        positionPointSizePopover();
        pointSizeSlider?.focus();
    };

    if (pointSizeBtn && pointSizePopover) {
        pointSizeBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            const isOpen = pointSizePopover.classList.contains('open');
            if (isOpen) {
                closePointSizePopover();
            } else {
                openPointSizePopover();
            }
        });

        document.addEventListener('click', (e) => {
            const target = e.target as Node | null;
            if (!target) return;
            if (pointSizePopover.contains(target)) return;
            if (pointSizeBtn.contains(target)) return;
            closePointSizePopover();
        });

        window.addEventListener('resize', () => {
            if (!pointSizePopover.classList.contains('open')) return;
            positionPointSizePopover();
        });

        document.addEventListener('keydown', (e) => {
            if (e.key !== 'Escape') return;
            closePointSizePopover();
        });
    }

    // Label font size popover (opened by long-press on T button)
    const setLabelFontSizeUi = (scale: number) => {
        if (!labelFontSizeValue) return;
        labelFontSizeValue.textContent = `${Math.round(scale * 100)}%`;
    };

    const applyLabelFontScale = (scale: number) => {
        map.setCategoryLabelFontScale(scale);
        const effective = map.getCategoryLabelFontScale();
        if (labelFontSizeSlider) {
            labelFontSizeSlider.value = String(effective);
        }
        setLabelFontSizeUi(effective);
    };

    if (labelFontSizeSlider) {
        const onInput = () => {
            const value = Number.parseFloat(labelFontSizeSlider!.value);
            if (!Number.isFinite(value)) return;
            applyLabelFontScale(value);
        };
        labelFontSizeSlider.addEventListener('input', onInput);
        onInput();
        labelFontSizeSlider.disabled = true;
    }

    const positionLabelFontSizePopover = () => {
        if (!labelFontSizePopover || !labelsBtn) return;
        const panel = labelsBtn.closest('.panel') as HTMLElement | null;
        if (!panel) return;

        const panelRect = panel.getBoundingClientRect();
        const btnRect = labelsBtn.getBoundingClientRect();

        const padding = 8;
        const gap = 8;

        const popoverWidth = labelFontSizePopover.offsetWidth;
        const popoverHeight = labelFontSizePopover.offsetHeight;

        let left = btnRect.left - panelRect.left;
        let top = btnRect.bottom - panelRect.top + gap;

        const maxLeft = panelRect.width - popoverWidth - padding;
        if (left > maxLeft) left = maxLeft;
        if (left < padding) left = padding;

        if (top + popoverHeight > panelRect.height - padding) {
            top = btnRect.top - panelRect.top - popoverHeight - gap;
        }
        if (top < padding) top = padding;

        labelFontSizePopover.style.left = `${left}px`;
        labelFontSizePopover.style.top = `${top}px`;
    };

    const closeLabelFontSizePopover = () => {
        if (!labelFontSizePopover || !labelsBtn) return;
        labelFontSizePopover.classList.remove('open');
        labelFontSizePopover.setAttribute('aria-hidden', 'true');
        if (labelFontSizeSlider) {
            labelFontSizeSlider.disabled = true;
        }
    };

    const openLabelFontSizePopover = () => {
        if (!labelFontSizePopover || !labelsBtn) return;
        labelLongPressJustHandled = true;
        labelFontSizePopover.classList.add('open');
        labelFontSizePopover.setAttribute('aria-hidden', 'false');
        if (labelFontSizeSlider) {
            labelFontSizeSlider.disabled = false;
        }
        positionLabelFontSizePopover();
        labelFontSizeSlider?.focus();
    };

    const LONG_PRESS_MS = 500;
    let labelLongPressTimer: ReturnType<typeof setTimeout> | null = null;

    const clearLabelLongPressTimer = () => {
        if (labelLongPressTimer !== null) {
            clearTimeout(labelLongPressTimer);
            labelLongPressTimer = null;
        }
    };

    if (labelsBtn && labelFontSizePopover) {
        labelsBtn.addEventListener('mousedown', (e) => {
            if (e.button !== 0) return;
            clearLabelLongPressTimer();
            labelLongPressTimer = setTimeout(() => {
                labelLongPressTimer = null;
                openLabelFontSizePopover();
            }, LONG_PRESS_MS);
        });
        labelsBtn.addEventListener('mouseup', clearLabelLongPressTimer);
        labelsBtn.addEventListener('mouseleave', clearLabelLongPressTimer);

        labelsBtn.addEventListener('touchstart', () => {
            clearLabelLongPressTimer();
            labelLongPressTimer = setTimeout(() => {
                labelLongPressTimer = null;
                openLabelFontSizePopover();
            }, LONG_PRESS_MS);
        }, { passive: true });
        labelsBtn.addEventListener('touchend', clearLabelLongPressTimer, { passive: true });
        labelsBtn.addEventListener('touchcancel', clearLabelLongPressTimer, { passive: true });

        document.addEventListener('click', (e) => {
            const target = e.target as Node | null;
            if (!target) return;
            if (labelFontSizePopover.contains(target)) return;
            if (labelsBtn.contains(target)) return;
            closeLabelFontSizePopover();
        });

        window.addEventListener('resize', () => {
            if (!labelFontSizePopover.classList.contains('open')) return;
            positionLabelFontSizePopover();
        });

        document.addEventListener('keydown', (e) => {
            if (e.key !== 'Escape') return;
            closeLabelFontSizePopover();
        });
    }

    if (resetBtn) {
        resetBtn.addEventListener('click', () => {
            map.resetView();
            if (pointSizeSlider) {
                pointSizeSlider.value = '1';
            }
            applyPointSize(1);
            closePointSizePopover();
            if (labelFontSizeSlider) {
                labelFontSizeSlider.value = '1';
            }
            applyLabelFontScale(1);
            closeLabelFontSizePopover();

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
            });

            onAfterReset?.();
        });
    }

    if (panLockBtn) {
        updatePanLockUi();
        panLockBtn.addEventListener('click', () => {
            map.setPanLocked(!map.isPanLocked());
            updatePanLockUi();
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

    if (deBtn) {
        deBtn.addEventListener('click', () => {
            const url = new URL(window.location.href);
            url.pathname = '/de.html';
            url.searchParams.set('dataset', datasetId);
            window.location.href = url.toString();
        });
    }

    const blastpBtn = document.getElementById('btn-blastp') as HTMLButtonElement | null;
    if (blastpBtn) {
        blastpBtn.addEventListener('click', () => {
            window.location.href = '/blastp.html';
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

// Mobile UI: floating action buttons + floating cards
function setupMobileUI(tabPanel: TabPanel, mapController: MapController): void {
    const mobileQuery = window.matchMedia('(max-width: 768px)');
    if (!mobileQuery.matches) return;

    const mapContainer = document.getElementById('map-container');
    const tabPanelContainer = document.getElementById('tab-panel-container');
    const toolsPanel = document.getElementById('tools-panel');
    if (!mapContainer || !tabPanelContainer || !toolsPanel) return;

    // Track which card is open: 'tabs' | 'tools' | null
    let openCard: 'tabs' | 'tools' | null = null;
    const leafletMap = mapController.getMap();

    // Create backdrop
    const backdrop = document.createElement('div');
    backdrop.className = 'mobile-backdrop';
    mapContainer.appendChild(backdrop);

    // Create FABs container
    const fabsContainer = document.createElement('div');
    fabsContainer.className = 'mobile-fabs';
    fabsContainer.innerHTML = `
        <button class="mobile-fab-btn" data-action="category" title="Category">
            <svg viewBox="0 0 24 24" width="20" height="20"><path fill="currentColor" d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-1 17.93c-3.95-.49-7-3.85-7-7.93 0-.62.08-1.21.21-1.79L9 15v1c0 1.1.9 2 2 2v1.93zm6.9-2.54c-.26-.81-1-1.39-1.9-1.39h-1v-3c0-.55-.45-1-1-1H8v-2h2c.55 0 1-.45 1-1V7h2c1.1 0 2-.9 2-2v-.41c2.93 1.19 5 4.06 5 7.41 0 2.08-.8 3.97-2.1 5.39z"/></svg>
        </button>
        <button class="mobile-fab-btn" data-action="expression" title="Expression">
            <svg viewBox="0 0 24 24" width="20" height="20"><path fill="currentColor" d="M7 18c-1.1 0-1.99.9-1.99 2S5.9 22 7 22s2-.9 2-2-.9-2-2-2zM1 2v2h2l3.6 7.59-1.35 2.45c-.16.28-.25.61-.25.96 0 1.1.9 2 2 2h12v-2H7.42c-.14 0-.25-.11-.25-.25l.03-.12.9-1.63h7.45c.75 0 1.41-.41 1.75-1.03l3.58-6.49c.08-.14.12-.31.12-.48 0-.55-.45-1-1-1H5.21l-.94-2H1zm16 16c-1.1 0-1.99.9-1.99 2s.89 2 1.99 2 2-.9 2-2-.9-2-2-2z"/></svg>
        </button>
        <button class="mobile-fab-btn" data-action="tools" title="Tools">
            <svg viewBox="0 0 24 24" width="20" height="20"><path fill="currentColor" d="M22.7 19l-9.1-9.1c.9-2.3.4-5-1.5-6.9-2-2-5-2.4-7.4-1.3L9 6 6 9 1.6 4.7C.4 7.1.9 10.1 2.9 12.1c1.9 1.9 4.6 2.4 6.9 1.5l9.1 9.1c.4.4 1 .4 1.4 0l2.3-2.3c.5-.4.5-1.1.1-1.4z"/></svg>
        </button>
    `;
    mapContainer.appendChild(fabsContainer);

    // Create tabs card
    const tabsCard = document.createElement('div');
    tabsCard.className = 'mobile-card mobile-card--tabs';
    tabsCard.innerHTML = `
        <div class="mobile-card__header">
            <span class="mobile-card__title"></span>
            <button class="mobile-card__close" aria-label="Close">&times;</button>
        </div>
        <div class="mobile-card__body"></div>
    `;
    mapContainer.appendChild(tabsCard);

    // Create tools card
    const toolsCard = document.createElement('div');
    toolsCard.className = 'mobile-card mobile-card--tools';
    toolsCard.innerHTML = `
        <div class="mobile-card__header">
            <span class="mobile-card__title">Tools</span>
            <button class="mobile-card__close" aria-label="Close">&times;</button>
        </div>
        <div class="mobile-card__body"></div>
    `;
    mapContainer.appendChild(toolsCard);

    // Reparent content into cards
    const tabsCardBody = tabsCard.querySelector('.mobile-card__body')!;
    const toolsCardBody = toolsCard.querySelector('.mobile-card__body')!;
    tabsCardBody.appendChild(tabPanelContainer);
    toolsCardBody.appendChild(toolsPanel);

    // Card title element for tabs card
    const tabsCardTitle = tabsCard.querySelector('.mobile-card__title') as HTMLElement;

    const updateTabsCardTitle = () => {
        const currentTab = tabPanel.getCurrentTab();
        tabsCardTitle.textContent = currentTab === 'category' ? 'Category' : 'Expression';
    };

    const closeAllCards = () => {
        tabsCard.classList.remove('open');
        toolsCard.classList.remove('open');
        backdrop.classList.remove('open');
        openCard = null;
        // Re-enable map dragging when cards are closed
        leafletMap.dragging.enable();
        leafletMap.touchZoom.enable();
    };

    const openTabsCard = (tab: 'category' | 'expression') => {
        tabPanel.switchTab(tab);
        updateTabsCardTitle();
        toolsCard.classList.remove('open');
        tabsCard.classList.add('open');
        backdrop.classList.add('open');
        openCard = 'tabs';
        // Disable map dragging to prevent touch events from moving the map
        leafletMap.dragging.disable();
        leafletMap.touchZoom.disable();
    };

    const openToolsCard = () => {
        tabsCard.classList.remove('open');
        toolsCard.classList.add('open');
        backdrop.classList.add('open');
        openCard = 'tools';
        // Disable map dragging to prevent touch events from moving the map
        leafletMap.dragging.disable();
        leafletMap.touchZoom.disable();
    };

    // FAB click handlers
    fabsContainer.addEventListener('click', (e) => {
        const btn = (e.target as HTMLElement).closest('.mobile-fab-btn') as HTMLElement | null;
        if (!btn) return;
        const action = btn.dataset.action;
        if (action === 'category') {
            if (openCard === 'tabs' && tabPanel.getCurrentTab() === 'category') {
                closeAllCards();
            } else {
                openTabsCard('category');
            }
        } else if (action === 'expression') {
            if (openCard === 'tabs' && tabPanel.getCurrentTab() === 'expression') {
                closeAllCards();
            } else {
                openTabsCard('expression');
            }
        } else if (action === 'tools') {
            if (openCard === 'tools') {
                closeAllCards();
            } else {
                openToolsCard();
            }
        }
    });

    // Close button handlers
    tabsCard.querySelector('.mobile-card__close')!.addEventListener('click', closeAllCards);
    toolsCard.querySelector('.mobile-card__close')!.addEventListener('click', closeAllCards);

    // Backdrop click closes cards
    backdrop.addEventListener('click', closeAllCards);

    // Escape key closes cards
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && openCard) {
            closeAllCards();
        }
    });

    // Hide .tab-header inside mobile card (tabs are controlled by FABs)
    const style = document.createElement('style');
    style.textContent = `
        .mobile-card--tabs .tab-header { display: none; }
        .mobile-card--tabs .tab-pane { display: block !important; padding: 0; }
        .mobile-card--tabs .tab-pane:not(.active) { display: none !important; }
    `;
    document.head.appendChild(style);
}

function formatNumber(n: number): string {
    if (n >= 1e6) return (n / 1e6).toFixed(1) + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(1) + 'K';
    return n.toString();
}

// Start application
document.addEventListener('DOMContentLoaded', init);
