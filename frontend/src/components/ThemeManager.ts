// Theme Manager - handles dark/light theme switching with localStorage persistence

export type Theme = 'dark' | 'light';

const STORAGE_KEY = 'soma-tiles-theme';

export class ThemeManager {
    private theme: Theme;

    constructor() {
        this.theme = this.loadTheme();
        this.applyTheme();
    }

    private loadTheme(): Theme {
        const stored = localStorage.getItem(STORAGE_KEY);
        return stored === 'light' ? 'light' : 'dark';
    }

    private applyTheme(): void {
        document.documentElement.dataset.theme = this.theme;
    }

    getTheme(): Theme {
        return this.theme;
    }

    setTheme(theme: Theme): void {
        this.theme = theme;
        localStorage.setItem(STORAGE_KEY, theme);
        this.applyTheme();
    }

    toggle(): Theme {
        const newTheme = this.theme === 'dark' ? 'light' : 'dark';
        this.setTheme(newTheme);
        return newTheme;
    }
}

