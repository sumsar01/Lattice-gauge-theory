/** Complex vectors are Float64Array of length 2n, interleaved (re, im). */
export type CVec = Float64Array;
export type MatVec = (x: CVec, y: CVec) => void;

export function popcount32(v: number): number {
  v = v - ((v >>> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return (((v + (v >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
}

export function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function norm(x: CVec): number {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i] * x[i];
  return Math.sqrt(s);
}

export function scale(x: CVec, a: number) {
  for (let i = 0; i < x.length; i++) x[i] *= a;
}

/** <x|y> */
export function dot(x: CVec, y: CVec): [number, number] {
  let re = 0;
  let im = 0;
  for (let i = 0; i < x.length; i += 2) {
    re += x[i] * y[i] + x[i + 1] * y[i + 1];
    im += x[i] * y[i + 1] - x[i + 1] * y[i];
  }
  return [re, im];
}

/** Real part of <x|y>, enough for Hermitian Rayleigh quotients. */
export function dotRe(x: CVec, y: CVec): number {
  let re = 0;
  for (let i = 0; i < x.length; i++) re += x[i] * y[i];
  return re;
}

/**
 * Eigen-decomposition of a symmetric tridiagonal matrix by implicit QL
 * (Numerical Recipes tqli). d: diagonal (overwritten with eigenvalues),
 * e: sub-diagonal with e[0] unused. Returns eigenvectors column-major-in-rows:
 * z[k * n + i] is component k of eigenvector i.
 */
export function tridiagEigen(dIn: ArrayLike<number>, eIn: ArrayLike<number>): { values: Float64Array; vectors: Float64Array } {
  const n = dIn.length;
  const d = Float64Array.from(dIn);
  const e = new Float64Array(n);
  for (let i = 1; i < n; i++) e[i - 1] = eIn[i];
  const z = new Float64Array(n * n);
  for (let i = 0; i < n; i++) z[i * n + i] = 1;

  for (let l = 0; l < n; l++) {
    let iter = 0;
    let m: number;
    do {
      for (m = l; m < n - 1; m++) {
        const dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
        if (Math.abs(e[m]) <= Number.EPSILON * dd) break;
      }
      if (m !== l) {
        if (iter++ === 100) throw new Error('tridiagEigen: no convergence');
        let g = (d[l + 1] - d[l]) / (2 * e[l]);
        let r = Math.hypot(g, 1);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? r : -r));
        let s = 1;
        let c = 1;
        let p = 0;
        let i: number;
        for (i = m - 1; i >= l; i--) {
          let f = s * e[i];
          const b = c * e[i];
          r = Math.hypot(f, g);
          e[i + 1] = r;
          if (r === 0) {
            d[i + 1] -= p;
            e[m] = 0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          for (let k = 0; k < n; k++) {
            f = z[k * n + i + 1];
            z[k * n + i + 1] = s * z[k * n + i] + c * f;
            z[k * n + i] = c * z[k * n + i] - s * f;
          }
        }
        if (r === 0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0;
      }
    } while (m !== l);
  }
  return { values: d, vectors: z };
}

function lowest(values: Float64Array): number {
  let k = 0;
  for (let i = 1; i < values.length; i++) if (values[i] < values[k]) k = i;
  return k;
}

/** One Lanczos run from v0; returns the Ritz ground state (two-pass, O(n) memory). */
function lanczosOnce(apply: MatVec, v0: CVec, maxIter: number, tol: number): { energy: number; state: CVec } {
  const n = v0.length;
  const alphas: number[] = [];
  const betas: number[] = [0];
  let vPrev = new Float64Array(n);
  let v = Float64Array.from(v0);
  scale(v, 1 / norm(v));
  let w = new Float64Array(n);
  let lastE = Infinity;

  for (let k = 0; k < maxIter; k++) {
    apply(v, w);
    const a = dotRe(v, w);
    alphas.push(a);
    const bPrev = betas[k];
    for (let i = 0; i < n; i++) w[i] -= a * v[i] + bPrev * vPrev[i];
    const b = norm(w);
    if ((k + 1) % 5 === 0 || b < 1e-12 || k === maxIter - 1) {
      const { values } = tridiagEigen(alphas, betas);
      const E = values[lowest(values)];
      if (Math.abs(E - lastE) < tol || b < 1e-12) break;
      lastE = E;
    }
    betas.push(b);
    const t = vPrev;
    vPrev = v;
    v = w;
    scale(v, 1 / b);
    w = t;
  }

  const m = alphas.length;
  const { values, vectors } = tridiagEigen(alphas, betas.slice(0, m));
  const k0 = lowest(values);

  const out = new Float64Array(n);
  vPrev.fill(0);
  v = Float64Array.from(v0);
  scale(v, 1 / norm(v));
  for (let k = 0; k < m; k++) {
    const c = vectors[k * m + k0];
    for (let i = 0; i < n; i++) out[i] += c * v[i];
    if (k === m - 1) break;
    apply(v, w);
    for (let i = 0; i < n; i++) w[i] -= alphas[k] * v[i] + betas[k] * vPrev[i];
    scale(w, 1 / betas[k + 1]);
    const t = vPrev;
    vPrev = v;
    v = w;
    w = t;
  }
  scale(out, 1 / norm(out));
  return { energy: values[k0], state: out };
}

export function groundState(
  apply: MatVec,
  dim: number,
  opts: { seed?: number; maxIter?: number; tol?: number; restarts?: number } = {},
): { energy: number; state: CVec; residual: number } {
  const rand = mulberry32(opts.seed ?? 12345);
  let v = new Float64Array(2 * dim);
  for (let i = 0; i < v.length; i++) v[i] = rand() - 0.5;
  const tmp = new Float64Array(2 * dim);
  let energy = 0;
  let residual = Infinity;
  for (let r = 0; r <= (opts.restarts ?? 4); r++) {
    const res = lanczosOnce(apply, v, Math.min(opts.maxIter ?? 200, 2 * dim), opts.tol ?? 1e-11);
    v = res.state;
    energy = res.energy;
    apply(v, tmp);
    for (let i = 0; i < tmp.length; i++) tmp[i] -= energy * v[i];
    residual = norm(tmp);
    if (residual < 1e-7) break;
  }
  return { energy, state: v, residual };
}

/**
 * x <- exp(-i H dt) x via a Krylov subspace of size m. Returns the
 * subspace error estimate (last beta times the last component).
 */
export function krylovStep(apply: MatVec, x: CVec, dt: number, m = 30): number {
  const n = x.length;
  const nx = norm(x);
  if (nx === 0) return 0;
  const V: Float64Array[] = [];
  const alphas: number[] = [];
  const betas: number[] = [0];
  let v = Float64Array.from(x);
  scale(v, 1 / nx);
  let w = new Float64Array(n);
  let lastBeta = 0;
  for (let k = 0; k < m; k++) {
    V.push(v);
    apply(v, w);
    const a = dotRe(v, w);
    alphas.push(a);
    for (let i = 0; i < n; i++) w[i] -= a * v[i];
    if (k > 0) {
      const vp = V[k - 1];
      for (let i = 0; i < n; i++) w[i] -= betas[k] * vp[i];
    }
    // one round of full reorthogonalisation keeps the small basis honest
    for (const u of V) {
      const [cr, ci] = dot(u, w);
      for (let i = 0; i < n; i += 2) {
        w[i] -= cr * u[i] - ci * u[i + 1];
        w[i + 1] -= cr * u[i + 1] + ci * u[i];
      }
    }
    const b = norm(w);
    lastBeta = b;
    if (b < 1e-12 || k === m - 1) break;
    betas.push(b);
    v = w;
    scale(v, 1 / b);
    w = new Float64Array(n);
  }

  const k = alphas.length;
  const { values, vectors } = tridiagEigen(alphas, betas.slice(0, k));
  // y = exp(-i T dt) e_0
  const yr = new Float64Array(k);
  const yi = new Float64Array(k);
  for (let j = 0; j < k; j++) {
    const c0 = vectors[0 * k + j];
    const cr = Math.cos(values[j] * dt) * c0;
    const ci = -Math.sin(values[j] * dt) * c0;
    for (let a = 0; a < k; a++) {
      yr[a] += vectors[a * k + j] * cr;
      yi[a] += vectors[a * k + j] * ci;
    }
  }
  x.fill(0);
  for (let a = 0; a < k; a++) {
    const u = V[a];
    const cr = yr[a] * nx;
    const ci = yi[a] * nx;
    for (let i = 0; i < n; i += 2) {
      x[i] += cr * u[i] - ci * u[i + 1];
      x[i + 1] += cr * u[i + 1] + ci * u[i];
    }
  }
  return lastBeta * Math.hypot(yr[k - 1], yi[k - 1]);
}
