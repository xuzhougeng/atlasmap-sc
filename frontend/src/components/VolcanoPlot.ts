export interface VolcanoPlotPoint {
    x: number; // log2fc
    y: number; // -log10(metric)
}

export interface VolcanoPlotConfig {
    pointSize?: number;
}

interface Rgb {
    r: number;
    g: number;
    b: number;
}

function clamp01(v: number): number {
    return Math.min(1, Math.max(0, v));
}

function parseCssColor(input: string): Rgb | null {
    const s = input.trim();
    if (!s) return null;

    // #rgb or #rrggbb
    if (s.startsWith('#')) {
        const hex = s.slice(1);
        if (hex.length === 3) {
            const r = parseInt(hex[0] + hex[0], 16);
            const g = parseInt(hex[1] + hex[1], 16);
            const b = parseInt(hex[2] + hex[2], 16);
            return { r, g, b };
        }
        if (hex.length === 6) {
            const r = parseInt(hex.slice(0, 2), 16);
            const g = parseInt(hex.slice(2, 4), 16);
            const b = parseInt(hex.slice(4, 6), 16);
            return { r, g, b };
        }
        return null;
    }

    // rgb(...) or rgba(...)
    const m = s.match(/^rgba?\((.+)\)$/i);
    if (m) {
        const parts = m[1].split(',').map(p => p.trim());
        if (parts.length < 3) return null;
        const r = Number(parts[0]);
        const g = Number(parts[1]);
        const b = Number(parts[2]);
        if (![r, g, b].every(Number.isFinite)) return null;
        return { r, g, b };
    }

    return null;
}

function rgbToVec3(c: Rgb): [number, number, number] {
    return [clamp01(c.r / 255), clamp01(c.g / 255), clamp01(c.b / 255)];
}

function createShader(gl: WebGLRenderingContext, type: number, source: string): WebGLShader {
    const shader = gl.createShader(type);
    if (!shader) {
        throw new Error('Failed to create WebGL shader');
    }
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    const ok = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
    if (!ok) {
        const log = gl.getShaderInfoLog(shader) || 'Unknown shader compile error';
        gl.deleteShader(shader);
        throw new Error(log);
    }
    return shader;
}

function createProgram(gl: WebGLRenderingContext, vsSource: string, fsSource: string): WebGLProgram {
    const vs = createShader(gl, gl.VERTEX_SHADER, vsSource);
    const fs = createShader(gl, gl.FRAGMENT_SHADER, fsSource);
    const program = gl.createProgram();
    if (!program) {
        gl.deleteShader(vs);
        gl.deleteShader(fs);
        throw new Error('Failed to create WebGL program');
    }
    gl.attachShader(program, vs);
    gl.attachShader(program, fs);
    gl.linkProgram(program);
    gl.deleteShader(vs);
    gl.deleteShader(fs);
    const ok = gl.getProgramParameter(program, gl.LINK_STATUS);
    if (!ok) {
        const log = gl.getProgramInfoLog(program) || 'Unknown program link error';
        gl.deleteProgram(program);
        throw new Error(log);
    }
    return program;
}

export class VolcanoPlot {
    private glCanvas: HTMLCanvasElement;
    private overlayCanvas: HTMLCanvasElement;
    private gl: WebGLRenderingContext;
    private overlay: CanvasRenderingContext2D;

    private program: WebGLProgram;
    private buffer: WebGLBuffer;
    private aXY: number;
    private uMin: WebGLUniformLocation;
    private uMax: WebGLUniformLocation;
    private uPointSize: WebGLUniformLocation;
    private uFcThreshold: WebGLUniformLocation;
    private uSigThreshold: WebGLUniformLocation;
    private uColorNeutral: WebGLUniformLocation;
    private uColorUp: WebGLUniformLocation;
    private uColorDown: WebGLUniformLocation;

    private points: Float32Array | null = null;
    private pointCount: number = 0;

    private cssWidth: number = 0;
    private cssHeight: number = 0;

    private minX: number = -1;
    private maxX: number = 1;
    private minY: number = 0;
    private maxY: number = 1;

    private xLabel: string = 'log2FC';
    private yLabel: string = '-log10(metric)';

    private pointSize: number;
    private fcThreshold: number = 1;
    private sigThreshold: number = 1.301;

    private resizeObserver: ResizeObserver | null = null;

    constructor(glCanvas: HTMLCanvasElement, overlayCanvas: HTMLCanvasElement, config: VolcanoPlotConfig = {}) {
        this.glCanvas = glCanvas;
        this.overlayCanvas = overlayCanvas;

        const gl = glCanvas.getContext('webgl', { antialias: true, alpha: true });
        if (!gl) {
            throw new Error('WebGL is not supported in this browser');
        }
        this.gl = gl;

        const overlay = overlayCanvas.getContext('2d');
        if (!overlay) {
            throw new Error('Failed to create 2D canvas context');
        }
        this.overlay = overlay;

        const vsSource = `
            attribute vec2 a_xy;

            uniform vec2 u_min;
            uniform vec2 u_max;
            uniform float u_pointSize;

            uniform float u_fcThreshold;
            uniform float u_sigThreshold;
            uniform vec3 u_colorNeutral;
            uniform vec3 u_colorUp;
            uniform vec3 u_colorDown;

            varying vec3 v_color;

            void main() {
                vec2 range = u_max - u_min;
                range = max(range, vec2(1e-6, 1e-6));
                vec2 norm = (a_xy - u_min) / range;   // 0..1
                vec2 clip = norm * 2.0 - 1.0;         // -1..1
                gl_Position = vec4(clip, 0.0, 1.0);
                gl_PointSize = u_pointSize;

                float x = a_xy.x;
                float y = a_xy.y;
                if (y >= u_sigThreshold && x >= u_fcThreshold) {
                    v_color = u_colorUp;
                } else if (y >= u_sigThreshold && x <= -u_fcThreshold) {
                    v_color = u_colorDown;
                } else {
                    v_color = u_colorNeutral;
                }
            }
        `;

        const fsSource = `
            precision mediump float;
            varying vec3 v_color;

            void main() {
                vec2 c = gl_PointCoord - vec2(0.5);
                float d = dot(c, c);
                if (d > 0.25) discard;
                gl_FragColor = vec4(v_color, 1.0);
            }
        `;

        this.program = createProgram(gl, vsSource, fsSource);

        const buffer = gl.createBuffer();
        if (!buffer) {
            gl.deleteProgram(this.program);
            throw new Error('Failed to create WebGL buffer');
        }
        this.buffer = buffer;

        this.aXY = gl.getAttribLocation(this.program, 'a_xy');
        const uMin = gl.getUniformLocation(this.program, 'u_min');
        const uMax = gl.getUniformLocation(this.program, 'u_max');
        const uPointSize = gl.getUniformLocation(this.program, 'u_pointSize');
        const uFcThreshold = gl.getUniformLocation(this.program, 'u_fcThreshold');
        const uSigThreshold = gl.getUniformLocation(this.program, 'u_sigThreshold');
        const uColorNeutral = gl.getUniformLocation(this.program, 'u_colorNeutral');
        const uColorUp = gl.getUniformLocation(this.program, 'u_colorUp');
        const uColorDown = gl.getUniformLocation(this.program, 'u_colorDown');
        if (!uMin || !uMax || !uPointSize || !uFcThreshold || !uSigThreshold || !uColorNeutral || !uColorUp || !uColorDown) {
            gl.deleteBuffer(this.buffer);
            gl.deleteProgram(this.program);
            throw new Error('Failed to locate WebGL uniforms');
        }
        this.uMin = uMin;
        this.uMax = uMax;
        this.uPointSize = uPointSize;
        this.uFcThreshold = uFcThreshold;
        this.uSigThreshold = uSigThreshold;
        this.uColorNeutral = uColorNeutral;
        this.uColorUp = uColorUp;
        this.uColorDown = uColorDown;

        this.pointSize = config.pointSize ?? 3;

        this.setupResizeObserver();
    }

    destroy(): void {
        if (this.resizeObserver) {
            this.resizeObserver.disconnect();
            this.resizeObserver = null;
        }
        this.gl.deleteBuffer(this.buffer);
        this.gl.deleteProgram(this.program);
    }

    setLabels(xLabel: string, yLabel: string): void {
        this.xLabel = xLabel;
        this.yLabel = yLabel;
        this.draw();
    }

    setThresholds(fcThreshold: number, sigThreshold: number): void {
        this.fcThreshold = Math.max(0, fcThreshold);
        this.sigThreshold = Math.max(0, sigThreshold);
        this.ensureRangesContainThresholds();
        this.draw();
    }

    setPointSize(size: number): void {
        this.pointSize = Math.max(1, size);
        this.draw();
    }

    setData(points: VolcanoPlotPoint[]): void {
        const xs: number[] = [];
        const ys: number[] = [];
        for (const p of points) {
            if (!Number.isFinite(p.x) || !Number.isFinite(p.y)) continue;
            xs.push(p.x);
            ys.push(p.y);
        }

        this.pointCount = xs.length;
        if (this.pointCount === 0) {
            this.points = null;
            this.minX = -1;
            this.maxX = 1;
            this.minY = 0;
            this.maxY = 1;
            this.draw();
            return;
        }

        const data = new Float32Array(this.pointCount * 2);
        let minX = Infinity;
        let maxX = -Infinity;
        let minY = Infinity;
        let maxY = -Infinity;
        for (let i = 0; i < this.pointCount; i++) {
            const x = xs[i];
            const y = ys[i];
            data[i * 2] = x;
            data[i * 2 + 1] = y;
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
        this.points = data;

        // Volcano plot conventions: symmetric x-range around 0, y starts at 0.
        const maxAbsX = Math.max(Math.abs(minX), Math.abs(maxX), this.fcThreshold * 1.2, 1e-6);
        this.minX = -maxAbsX;
        this.maxX = maxAbsX;
        this.minY = 0;
        this.maxY = Math.max(maxY, this.sigThreshold * 1.2, 1);

        this.upload();
        this.draw();
    }

    private ensureRangesContainThresholds(): void {
        const maxAbsX = Math.max(Math.abs(this.minX), Math.abs(this.maxX), this.fcThreshold * 1.2, 1e-6);
        this.minX = -maxAbsX;
        this.maxX = maxAbsX;
        this.maxY = Math.max(this.maxY, this.sigThreshold * 1.2, 1);
    }

    private upload(): void {
        if (!this.points) return;
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.buffer);
        this.gl.bufferData(this.gl.ARRAY_BUFFER, this.points, this.gl.STATIC_DRAW);
    }

    private setupResizeObserver(): void {
        const parent = this.glCanvas.parentElement;
        if (!parent || typeof ResizeObserver === 'undefined') return;

        this.resizeObserver = new ResizeObserver(() => {
            this.draw();
        });
        this.resizeObserver.observe(parent);
    }

    private resizeToDisplaySize(): boolean {
        const rect = this.glCanvas.getBoundingClientRect();
        const cssW = Math.floor(rect.width);
        const cssH = Math.floor(rect.height);
        if (cssW <= 0 || cssH <= 0) return false;

        const dpr = Math.max(1, window.devicePixelRatio || 1);
        const w = Math.max(1, Math.floor(cssW * dpr));
        const h = Math.max(1, Math.floor(cssH * dpr));

        let changed = false;
        if (this.glCanvas.width !== w || this.glCanvas.height !== h) {
            this.glCanvas.width = w;
            this.glCanvas.height = h;
            changed = true;
        }
        if (this.overlayCanvas.width !== w || this.overlayCanvas.height !== h) {
            this.overlayCanvas.width = w;
            this.overlayCanvas.height = h;
            changed = true;
        }

        this.cssWidth = cssW;
        this.cssHeight = cssH;
        return changed;
    }

    private getThemeColors(): { neutral: [number, number, number]; up: [number, number, number]; down: [number, number, number] } {
        const styles = getComputedStyle(document.documentElement);
        const neutralRaw = styles.getPropertyValue('--text-secondary');
        const upRaw = styles.getPropertyValue('--accent-primary');
        const downRaw = styles.getPropertyValue('--accent-secondary');
        const neutral = parseCssColor(neutralRaw) ?? { r: 160, g: 160, b: 160 };
        const up = parseCssColor(upRaw) ?? { r: 233, g: 69, b: 96 };
        const down = parseCssColor(downRaw) ?? { r: 15, g: 76, b: 117 };
        return { neutral: rgbToVec3(neutral), up: rgbToVec3(up), down: rgbToVec3(down) };
    }

    private drawOverlay(): void {
        const ctx = this.overlay;
        const dpr = Math.max(1, window.devicePixelRatio || 1);
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, this.cssWidth, this.cssHeight);

        const styles = getComputedStyle(document.documentElement);
        const borderColor = styles.getPropertyValue('--border-color').trim() || '#2a2a4a';
        const textSecondary = styles.getPropertyValue('--text-secondary').trim() || '#a0a0a0';

        const padL = 48;
        const padR = 16;
        const padT = 16;
        const padB = 36;

        const left = padL;
        const right = this.cssWidth - padR;
        const top = padT;
        const bottom = this.cssHeight - padB;

        const width = Math.max(1, right - left);
        const height = Math.max(1, bottom - top);

        const xToPx = (x: number) => {
            const t = (x - this.minX) / (this.maxX - this.minX);
            return left + t * width;
        };
        const yToPx = (y: number) => {
            const t = (y - this.minY) / (this.maxY - this.minY);
            return bottom - t * height;
        };

        ctx.save();
        ctx.strokeStyle = borderColor;
        ctx.lineWidth = 1;
        ctx.strokeRect(left, top, width, height);

        // Axes at x=0 and y=0 when within range.
        ctx.setLineDash([4, 4]);
        if (this.minX <= 0 && this.maxX >= 0) {
            const x0 = xToPx(0);
            ctx.beginPath();
            ctx.moveTo(x0, top);
            ctx.lineTo(x0, bottom);
            ctx.stroke();
        }
        if (this.minY <= 0 && this.maxY >= 0) {
            const y0 = yToPx(0);
            ctx.beginPath();
            ctx.moveTo(left, y0);
            ctx.lineTo(right, y0);
            ctx.stroke();
        }

        // Threshold lines.
        ctx.setLineDash([6, 4]);
        const xPos = xToPx(this.fcThreshold);
        const xNeg = xToPx(-this.fcThreshold);
        ctx.beginPath();
        ctx.moveTo(xPos, top);
        ctx.lineTo(xPos, bottom);
        ctx.moveTo(xNeg, top);
        ctx.lineTo(xNeg, bottom);
        ctx.stroke();

        const ySig = yToPx(this.sigThreshold);
        ctx.beginPath();
        ctx.moveTo(left, ySig);
        ctx.lineTo(right, ySig);
        ctx.stroke();
        ctx.setLineDash([]);

        // Labels
        ctx.fillStyle = textSecondary;
        ctx.font = '12px ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(this.xLabel, left + width / 2, this.cssHeight - 12);

        ctx.textAlign = 'left';
        ctx.fillText(this.yLabel, 8, 14);

        // Legend
        const colors = this.getThemeColors();
        const legendX = left + 8;
        const legendY = top + 10;
        const itemH = 14;
        const drawLegendItem = (idx: number, color: [number, number, number], label: string) => {
            const y = legendY + idx * itemH;
            ctx.fillStyle = `rgb(${Math.round(color[0] * 255)}, ${Math.round(color[1] * 255)}, ${Math.round(color[2] * 255)})`;
            ctx.beginPath();
            ctx.arc(legendX, y, 4, 0, Math.PI * 2);
            ctx.fill();
            ctx.fillStyle = textSecondary;
            ctx.fillText(label, legendX + 10, y + 4);
        };
        ctx.textAlign = 'left';
        drawLegendItem(0, colors.up, 'up');
        drawLegendItem(1, colors.down, 'down');
        drawLegendItem(2, colors.neutral, 'n.s.');

        ctx.restore();
    }

    draw(): void {
        const hasSize = this.resizeToDisplaySize();
        if (!hasSize && (this.glCanvas.width === 0 || this.glCanvas.height === 0)) {
            return;
        }

        const gl = this.gl;
        gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
        gl.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        if (!this.points || this.pointCount === 0) {
            this.drawOverlay();
            return;
        }

        gl.useProgram(this.program);

        gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);
        gl.enableVertexAttribArray(this.aXY);
        gl.vertexAttribPointer(this.aXY, 2, gl.FLOAT, false, 0, 0);

        gl.uniform2f(this.uMin, this.minX, this.minY);
        gl.uniform2f(this.uMax, this.maxX, this.maxY);

        const dpr = Math.max(1, window.devicePixelRatio || 1);
        gl.uniform1f(this.uPointSize, this.pointSize * dpr);
        gl.uniform1f(this.uFcThreshold, this.fcThreshold);
        gl.uniform1f(this.uSigThreshold, this.sigThreshold);

        const colors = this.getThemeColors();
        gl.uniform3fv(this.uColorNeutral, colors.neutral);
        gl.uniform3fv(this.uColorUp, colors.up);
        gl.uniform3fv(this.uColorDown, colors.down);

        gl.drawArrays(gl.POINTS, 0, this.pointCount);

        this.drawOverlay();
    }
}

