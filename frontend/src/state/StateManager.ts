// State Manager for SOMA-Tiles

export interface AppState {
    // View state
    zoom: number;
    center: [number, number];

    // Filter state
    selectedCategories: string[];
    expressionFilter: {
        gene: string | null;
        min: number;
        max: number;
    };

    // Selection state
    selection: {
        id: string | null;
        type: 'lasso' | 'rect' | null;
        polygon: [number, number][] | null;
        cellCount: number;
    };

    // Display state
    colorMode: 'category' | 'expression';
    colorGene: string | null;
    colorScale: string;
}

type StateListener = (state: AppState, prevState: AppState) => void;

export class StateManager {
    private state: AppState;
    private listeners: Map<string, Set<StateListener>>;

    constructor(initialState: AppState) {
        this.state = { ...initialState };
        this.listeners = new Map();
    }

    /**
     * Get current state (read-only)
     */
    getState(): Readonly<AppState> {
        return this.state;
    }

    /**
     * Update state with partial values
     */
    setState(partial: Partial<AppState>): void {
        const prevState = { ...this.state };
        this.state = { ...this.state, ...partial };
        this.notifyListeners(prevState);
    }

    /**
     * Subscribe to state changes
     */
    subscribe(key: string, callback: StateListener): () => void {
        if (!this.listeners.has(key)) {
            this.listeners.set(key, new Set());
        }
        this.listeners.get(key)!.add(callback);

        // Return unsubscribe function
        return () => {
            this.listeners.get(key)?.delete(callback);
        };
    }

    /**
     * Subscribe to all state changes
     */
    subscribeAll(callback: StateListener): () => void {
        return this.subscribe('*', callback);
    }

    private notifyListeners(prevState: AppState): void {
        // Notify specific listeners based on what changed
        const changedKeys = this.getChangedKeys(prevState, this.state);

        for (const key of changedKeys) {
            const keyListeners = this.listeners.get(key);
            if (keyListeners) {
                keyListeners.forEach(cb => cb(this.state, prevState));
            }
        }

        // Notify global listeners
        const globalListeners = this.listeners.get('*');
        if (globalListeners) {
            globalListeners.forEach(cb => cb(this.state, prevState));
        }
    }

    private getChangedKeys(prev: AppState, next: AppState): string[] {
        const keys: string[] = [];

        if (prev.zoom !== next.zoom) keys.push('zoom');
        if (prev.center !== next.center) keys.push('center');
        if (prev.selectedCategories !== next.selectedCategories) keys.push('selectedCategories');
        if (prev.expressionFilter !== next.expressionFilter) keys.push('expressionFilter');
        if (prev.selection !== next.selection) keys.push('selection');
        if (prev.colorMode !== next.colorMode) keys.push('colorMode');
        if (prev.colorGene !== next.colorGene) keys.push('colorGene');
        if (prev.colorScale !== next.colorScale) keys.push('colorScale');

        return keys;
    }
}
