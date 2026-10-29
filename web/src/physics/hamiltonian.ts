import { type Basis, linkCode } from './basis';
import { popcount32 } from './linalg';

/**
 * Kogut-Susskind Hamiltonian with two-component Wilson fermions (a = 1):
 *
 *   H = M sum_s (n_up - n_down)                      M = m + d r
 *     + g^2/2 sum_l E_l^2
 *     - 1/(2 g^2) sum_p (U_p + U_p^dag)
 *     - 1/2 sum_{s,k} psi^dag(s) (r sigma_z - i alpha_k) U_k(s) psi(s+k) + h.c.
 *
 * with alpha_x = sigma_x, alpha_y = sigma_y and U the |E| <= 1 truncated
 * raising operator. Fermions get Jordan-Wigner signs in mode order 2s + c.
 */
export interface Params {
  mass: number;
  g: number;
  r: number;
}

export interface CSR {
  rowPtr: Int32Array;
  col: Int32Array;
  re: Float64Array;
  im: Float64Array;
  /** plaquette or link id each entry belongs to */
  tag: Int32Array;
}

export interface HamiltonianParts {
  basis: Basis;
  dim: number;
  dims: 1 | 2;
  /** sum_s (n_up - n_down) */
  diagMass: Float64Array;
  /** sum_l E_l^2 */
  diagElectric: Float64Array;
  /** sum_p (U_p + U_p^dag), tag = plaquette */
  plaq: CSR;
  /** r-independent hopping, tag = link */
  hopA: CSR;
  /** hopping proportional to r, tag = link */
  hopR: CSR;
}

class CsrBuilder {
  private n = 0;
  private col = new Int32Array(1024);
  private re = new Float64Array(1024);
  private im = new Float64Array(1024);
  private tag = new Int32Array(1024);
  private rows: number[] = [0];

  push(j: number, re: number, im: number, tag: number) {
    if (this.n === this.col.length) this.grow();
    this.col[this.n] = j;
    this.re[this.n] = re;
    this.im[this.n] = im;
    this.tag[this.n] = tag;
    this.n++;
  }

  endRow() {
    this.rows.push(this.n);
  }

  private grow() {
    const cap = this.col.length * 2;
    const g = <T extends Int32Array | Float64Array>(a: T, C: new (n: number) => T): T => {
      const b = new C(cap);
      b.set(a);
      return b;
    };
    this.col = g(this.col, Int32Array);
    this.re = g(this.re, Float64Array);
    this.im = g(this.im, Float64Array);
    this.tag = g(this.tag, Int32Array);
  }

  build(): CSR {
    return {
      rowPtr: Int32Array.from(this.rows),
      col: this.col.slice(0, this.n),
      re: this.re.slice(0, this.n),
      im: this.im.slice(0, this.n),
      tag: this.tag.slice(0, this.n),
    };
  }
}

// [a, b, re, im] of -1/2 A_ab for the r-independent part, per direction
const HOP_A: number[][][] = [
  [
    [0, 1, 0, 0.5],
    [1, 0, 0, 0.5],
  ],
  [
    [0, 1, 0.5, 0],
    [1, 0, -0.5, 0],
  ],
];
const HOP_R: number[][] = [
  [0, 0, -0.5, 0],
  [1, 1, 0.5, 0],
];

const POW3: number[] = [];
for (let i = 0, p = 1; i <= 33; i++, p *= 3) POW3.push(p);

export function buildHamiltonian(basis: Basis): HamiltonianParts {
  const L = basis.lattice;
  const { nSites, nLinks, linkFrom, linkTo, linkDir, plaquettes, nPlaquettes } = L;
  const { dim, links, ferm, linkCodeStart } = basis;

  const lookup = (code: number, f: number): number => {
    const start = linkCodeStart.get(code);
    if (start === undefined) return -1;
    let rank = 0;
    let k = 0;
    for (let s = 0; s < nSites; s++) {
      const pair = (f >>> (2 * s)) & 3;
      if (pair === 1) rank |= 1 << k;
      if (pair === 1 || pair === 2) k++;
    }
    const idx = start + rank;
    return ferm[idx] === f ? idx : -1;
  };
  const mustFind = (code: number, f: number): number => {
    const j = lookup(code, f);
    if (j < 0) throw new Error('Gauss law violated by Hamiltonian term');
    return j;
  };

  const diagMass = new Float64Array(dim);
  const diagElectric = new Float64Array(dim);
  const plaq = new CsrBuilder();
  const hopA = new CsrBuilder();
  const hopR = new CsrBuilder();

  for (let i = 0; i < dim; i++) {
    const off = i * nLinks;
    const f = ferm[i];
    const code = linkCode(links, off, nLinks);

    let m = 0;
    for (let s = 0; s < nSites; s++) m += ((f >>> (2 * s)) & 1) - ((f >>> (2 * s + 1)) & 1);
    diagMass[i] = m;
    let e2 = 0;
    for (let l = 0; l < nLinks; l++) e2 += links[off + l] * links[off + l];
    diagElectric[i] = e2;

    for (let p = 0; p < nPlaquettes; p++) {
      const b = plaquettes[4 * p];
      const r = plaquettes[4 * p + 1];
      const t = plaquettes[4 * p + 2];
      const lft = plaquettes[4 * p + 3];
      const Eb = links[off + b];
      const Er = links[off + r];
      const Et = links[off + t];
      const El = links[off + lft];
      if (Eb < 1 && Er < 1 && Et > -1 && El > -1) {
        plaq.push(mustFind(code + POW3[b] + POW3[r] - POW3[t] - POW3[lft], f), 1, 0, p);
      }
      if (Eb > -1 && Er > -1 && Et < 1 && El < 1) {
        plaq.push(mustFind(code - POW3[b] - POW3[r] + POW3[t] + POW3[lft], f), 1, 0, p);
      }
    }
    plaq.endRow();

    for (let l = 0; l < nLinks; l++) {
      const El = links[off + l];
      const s = linkFrom[l];
      const t = linkTo[l];
      const d = linkDir[l];
      const hop = (b: CsrBuilder, a: number, bb: number, cre: number, cim: number) => {
        const mu = 2 * s + a;
        const nu = 2 * t + bb;
        const nMu = (f >>> mu) & 1;
        const nNu = (f >>> nu) & 1;
        if (nMu === nNu) return;
        const lo = Math.min(mu, nu);
        const hi = Math.max(mu, nu);
        const between = hi - lo - 1;
        const sign = between > 0 && popcount32((f >>> (lo + 1)) & ((1 << between) - 1)) & 1 ? -1 : 1;
        const f2 = (f ^ (1 << mu) ^ (1 << nu)) >>> 0;
        if (nNu === 1 && El < 1) {
          // c psi_mu^dag U psi_nu: H_ji = c sign, so row i holds conj(c) sign
          b.push(mustFind(code + POW3[l], f2), sign * cre, -sign * cim, l);
        } else if (nMu === 1 && El > -1) {
          // conj(c) psi_nu^dag U^dag psi_mu
          b.push(mustFind(code - POW3[l], f2), sign * cre, sign * cim, l);
        }
      };
      for (const [a, bb, cre, cim] of HOP_A[d]) hop(hopA, a, bb, cre, cim);
      for (const [a, bb, cre, cim] of HOP_R) hop(hopR, a, bb, cre, cim);
    }
    hopA.endRow();
    hopR.endRow();
  }

  return {
    basis,
    dim,
    dims: L.Ly > 1 ? 2 : 1,
    diagMass,
    diagElectric,
    plaq: plaq.build(),
    hopA: hopA.build(),
    hopR: hopR.build(),
  };
}

export const massGap = (H: HamiltonianParts, p: Params) => p.mass + H.dims * p.r;

function spmvAdd(A: CSR, w: number, x: Float64Array, y: Float64Array) {
  const { rowPtr, col, re, im } = A;
  const n = rowPtr.length - 1;
  for (let i = 0; i < n; i++) {
    let sr = 0;
    let si = 0;
    for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
      const j = col[k];
      const xr = x[2 * j];
      const xi = x[2 * j + 1];
      sr += re[k] * xr - im[k] * xi;
      si += re[k] * xi + im[k] * xr;
    }
    y[2 * i] += w * sr;
    y[2 * i + 1] += w * si;
  }
}

/** y = H x for complex vectors stored interleaved (re, im). */
export function applyH(H: HamiltonianParts, p: Params, x: Float64Array, y: Float64Array) {
  const M = massGap(H, p);
  const eW = (p.g * p.g) / 2;
  const n = H.dim;
  for (let i = 0; i < n; i++) {
    const d = M * H.diagMass[i] + eW * H.diagElectric[i];
    y[2 * i] = d * x[2 * i];
    y[2 * i + 1] = d * x[2 * i + 1];
  }
  spmvAdd(H.plaq, -1 / (2 * p.g * p.g), x, y);
  spmvAdd(H.hopA, 1, x, y);
  if (p.r !== 0) spmvAdd(H.hopR, p.r, x, y);
}

/** Re <x|A|x> split by entry tag. */
export function taggedExpectation(A: CSR, x: Float64Array, nTags: number): Float64Array {
  const out = new Float64Array(nTags);
  const { rowPtr, col, re, im, tag } = A;
  const n = rowPtr.length - 1;
  for (let i = 0; i < n; i++) {
    const ar = x[2 * i];
    const ai = x[2 * i + 1];
    if (ar === 0 && ai === 0) continue;
    for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
      const j = col[k];
      const br = x[2 * j];
      const bi = x[2 * j + 1];
      // Re[conj(a) (re + i im) b]
      const vr = re[k] * br - im[k] * bi;
      const vi = re[k] * bi + im[k] * br;
      out[tag[k]] += ar * vr + ai * vi;
    }
  }
  return out;
}
