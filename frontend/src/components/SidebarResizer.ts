// Sidebar Resizer Component
// Allows users to resize the sidebar by dragging

export class SidebarResizer {
    private resizer: HTMLElement;
    private sidebar: HTMLElement;
    private minWidth = 200;
    private maxWidth = 500;
    private storageKey = 'soma-tiles-sidebar-width';

    constructor(resizer: HTMLElement, sidebar: HTMLElement) {
        this.resizer = resizer;
        this.sidebar = sidebar;
        this.setupEvents();
        this.loadSavedWidth();
    }

    private setupEvents(): void {
        this.resizer.addEventListener('mousedown', this.startResize);
    }

    private startResize = (e: MouseEvent): void => {
        e.preventDefault();
        this.resizer.classList.add('active');
        document.body.style.cursor = 'col-resize';
        document.body.style.userSelect = 'none';
        document.addEventListener('mousemove', this.resize);
        document.addEventListener('mouseup', this.stopResize);
    };

    private resize = (e: MouseEvent): void => {
        const newWidth = Math.min(this.maxWidth, Math.max(this.minWidth, e.clientX));
        this.sidebar.style.width = `${newWidth}px`;
    };

    private stopResize = (): void => {
        this.resizer.classList.remove('active');
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
        document.removeEventListener('mousemove', this.resize);
        document.removeEventListener('mouseup', this.stopResize);
        this.saveWidth();
    };

    private saveWidth(): void {
        const width = this.sidebar.style.width;
        if (width) {
            localStorage.setItem(this.storageKey, width);
        }
    }

    private loadSavedWidth(): void {
        const saved = localStorage.getItem(this.storageKey);
        if (saved) {
            const width = parseInt(saved, 10);
            if (width >= this.minWidth && width <= this.maxWidth) {
                this.sidebar.style.width = saved;
            }
        }
    }
}
