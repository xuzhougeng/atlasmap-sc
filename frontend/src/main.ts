// SOMA-Tiles Frontend Entry Point

// Import CSS (Leaflet must come first so app styles can override)
import 'leaflet/dist/leaflet.css';
import './styles/main.css';

import { MapController } from './map/MapController';
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
    console.log('Initializing SOMA-Tiles...');
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
    const maxNativeZoom = Math.max(0, metadata.zoom_levels - 1);
    // Allow a few extra zoom levels to "magnify" the highest native tiles.
    const maxZoom = maxNativeZoom + 4;
    const mapController = new MapController(mapContainer, {
        apiUrl: tilesBaseUrl,
        coord: effectiveCoord ?? undefined,
        tileSize: 256,
        maxZoom,
        maxNativeZoom,
        // Start at zoom level 2 for better initial detail.
        initialZoom: 2,
        bounds: metadata.bounds,
    });

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
            showCategoryLabels = !showCategoryLabels;
            labelsBtn.classList.toggle('active', showCategoryLabels);
            void refreshCategoryLabels();
        });
    }

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
                } else {
                    categoryLegend.hide();
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
        const statsZoom = maxNativeZoom;
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
    cellQueryPanel.setMaxNativeZoom(maxNativeZoom);
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
    }

    // Initialize with default category (skip if gene was specified in URL)
    if (defaultCategory && !geneFromUrl) {
        const filter = categoryFilter?.getFilterForColumn(defaultCategory) ?? null;
        mapController.setCategoryColumn(defaultCategory, filter);
        categoryLegend.loadLegend(defaultCategory);
        currentCategoryFilter = filter ?? null;
        void refreshCategoryLabels();
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

    if (resetBtn) {
        resetBtn.addEventListener('click', () => {
            map.resetView();
            if (pointSizeSlider) {
                pointSizeSlider.value = '1';
            }
            applyPointSize(1);
            closePointSizePopover();

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

function formatNumber(n: number): string {
    if (n >= 1e6) return (n / 1e6).toFixed(1) + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(1) + 'K';
    return n.toString();
}

// Start application
document.addEventListener('DOMContentLoaded', init);
