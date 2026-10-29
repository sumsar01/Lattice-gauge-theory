import { SERIES } from '../viz/colors';

interface Panel {
  title: string;
  color: string;
  values: number[];
}

const INK = { primary: '#e8eaf6', secondary: '#9aa0bd', muted: '#5d6384', grid: 'rgba(154,160,189,0.14)' };

function niceCeil(v: number): number {
  const p = 10 ** Math.floor(Math.log10(v));
  for (const m of [1, 2, 4, 5, 10]) if (m * p >= v) return m * p;
  return 10 * p;
}

/** Two stacked small-multiple line charts sharing one time axis. */
export class TimeChart {
  private canvas = document.createElement('canvas');
  private tooltip = document.createElement('div');
  private times: number[] = [];
  private tMax = 1;
  private playhead = 0;
  private hoverX: number | null = null;
  private panels: Panel[] = [
    { title: 'Particles + antiparticles', color: SERIES.particles, values: [] },
    { title: 'Electric field energy', color: SERIES.field, values: [] },
  ];
  onSeek: ((t: number) => void) | null = null;

  constructor(private host: HTMLElement) {
    host.classList.add('chart');
    this.tooltip.className = 'chart-tip';
    host.append(this.canvas, this.tooltip);
    new ResizeObserver(() => this.draw()).observe(host);
    this.canvas.addEventListener('pointermove', (e) => {
      this.hoverX = e.offsetX;
      this.draw();
    });
    this.canvas.addEventListener('pointerleave', () => {
      this.hoverX = null;
      this.draw();
    });
    this.canvas.addEventListener('click', (e) => {
      const t = this.timeAt(e.offsetX);
      if (t !== null) this.onSeek?.(t);
    });
  }

  reset(tMax: number) {
    this.tMax = tMax;
    this.times = [];
    for (const p of this.panels) p.values = [];
    this.draw();
  }

  push(t: number, particles: number, field: number) {
    this.times.push(t);
    this.panels[0].values.push(particles);
    this.panels[1].values.push(field);
  }

  setPlayhead(t: number) {
    this.playhead = t;
    this.draw();
  }

  private readonly pad = { l: 34, r: 10, t: 22, b: 20, gap: 26 };

  private timeAt(x: number): number | null {
    const w = this.host.clientWidth;
    const { l, r } = this.pad;
    if (x < l || x > w - r) return null;
    return ((x - l) / (w - l - r)) * this.tMax;
  }

  private indexAt(t: number): number {
    if (!this.times.length) return -1;
    let lo = 0;
    let hi = this.times.length - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (this.times[mid] <= t) lo = mid;
      else hi = mid - 1;
    }
    return lo;
  }

  draw() {
    const w = this.host.clientWidth;
    const h = this.host.clientHeight;
    if (!w || !h) return;
    const dpr = Math.min(window.devicePixelRatio, 2);
    this.canvas.width = w * dpr;
    this.canvas.height = h * dpr;
    this.canvas.style.width = `${w}px`;
    this.canvas.style.height = `${h}px`;
    const g = this.canvas.getContext('2d')!;
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);

    const { l, r, t, b, gap } = this.pad;
    const ph = (h - t - gap - b) / 2;
    const xOf = (time: number) => l + (time / this.tMax) * (w - l - r);
    const hoverT = this.hoverX !== null ? this.timeAt(this.hoverX) : null;
    const hoverIdx = hoverT !== null ? this.indexAt(hoverT) : -1;

    this.panels.forEach((p, k) => {
      const top = t + k * (ph + gap);
      const max = niceCeil(Math.max(1, ...p.values));
      const yOf = (v: number) => top + ph - (v / max) * ph;

      g.font = '600 11px Inter, system-ui, sans-serif';
      g.fillStyle = INK.secondary;
      g.textBaseline = 'bottom';
      g.fillText(p.title, l, top - 6);

      g.strokeStyle = INK.grid;
      g.lineWidth = 1;
      g.fillStyle = INK.muted;
      g.font = '10px Inter, system-ui, sans-serif';
      g.textBaseline = 'middle';
      g.textAlign = 'right';
      for (const frac of [0, 0.5, 1]) {
        const v = max * frac;
        const y = yOf(v);
        g.beginPath();
        g.moveTo(l, y);
        g.lineTo(w - r, y);
        g.stroke();
        g.fillText(String(+v.toFixed(1)), l - 6, y);
      }
      g.textAlign = 'left';

      if (p.values.length > 1) {
        g.strokeStyle = p.color;
        g.lineWidth = 2;
        g.lineJoin = 'round';
        g.beginPath();
        p.values.forEach((v, i) => {
          const x = xOf(this.times[i]);
          if (i === 0) g.moveTo(x, yOf(v));
          else g.lineTo(x, yOf(v));
        });
        g.stroke();
      }

      if (hoverIdx >= 0) {
        const x = xOf(this.times[hoverIdx]);
        const y = yOf(p.values[hoverIdx]);
        g.fillStyle = p.color;
        g.strokeStyle = SERIES.surface;
        g.lineWidth = 2;
        g.beginPath();
        g.arc(x, y, 4.5, 0, Math.PI * 2);
        g.fill();
        g.stroke();
      }
    });

    g.fillStyle = INK.muted;
    g.font = '10px Inter, system-ui, sans-serif';
    g.textBaseline = 'alphabetic';
    for (let k = 0; k <= 4; k++) {
      const time = (this.tMax * k) / 4;
      g.textAlign = k === 0 ? 'left' : k === 4 ? 'right' : 'center';
      g.fillText(k === 4 ? `${time.toFixed(0)} time` : time.toFixed(0), xOf(time), h - 4);
    }
    g.textAlign = 'left';

    const px = xOf(Math.min(this.playhead, this.tMax));
    g.strokeStyle = 'rgba(232,234,246,0.55)';
    g.lineWidth = 1;
    g.beginPath();
    g.moveTo(px, t - 4);
    g.lineTo(px, h - b);
    g.stroke();

    if (hoverIdx >= 0 && this.hoverX !== null) {
      const x = xOf(this.times[hoverIdx]);
      g.strokeStyle = 'rgba(232,234,246,0.25)';
      g.setLineDash([3, 3]);
      g.beginPath();
      g.moveTo(x, t - 4);
      g.lineTo(x, h - b);
      g.stroke();
      g.setLineDash([]);
      const [a, e] = this.panels.map((p) => p.values[hoverIdx]);
      this.tooltip.innerHTML =
        `<div class="tip-t">time ${this.times[hoverIdx].toFixed(2)}</div>` +
        `<div><i style="background:${SERIES.particles}"></i>particles <b>${a.toFixed(2)}</b></div>` +
        `<div><i style="background:${SERIES.field}"></i>field energy <b>${e.toFixed(2)}</b></div>` +
        `<div class="tip-hint">click to jump here</div>`;
      this.tooltip.style.display = 'block';
      const tw = this.tooltip.offsetWidth;
      this.tooltip.style.left = `${Math.min(w - tw - 4, Math.max(4, x + 12))}px`;
      this.tooltip.style.top = `4px`;
    } else {
      this.tooltip.style.display = 'none';
    }
  }
}
