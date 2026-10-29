import type { Lattice } from './lattice';

/**
 * Gauge-invariant basis of the truncated (|E| <= 1) U(1) lattice gauge theory.
 *
 * Each site s carries two fermion modes (Dirac components), mode index 2s + c.
 * Dynamical charge q(s) = n(2s) + n(2s+1) - 1, so the Dirac-sea vacuum is
 * n(2s) = 0, n(2s+1) = 1. Gauss law: sum_out E - sum_in E = q(s) + Q(s),
 * with Q a static background charge.
 *
 * Instead of filtering the full Hilbert space we enumerate link
 * configurations (pruned depth-first) and solve Gauss law for the fermions:
 * q = +-1 fixes both modes, q = 0 leaves a two-fold choice.
 */
export interface Basis {
  lattice: Lattice;
  staticCharges: Int8Array;
  dim: number;
  /** E on link l of state i at [i * nLinks + l], values in {-1, 0, 1} */
  links: Int8Array;
  /** fermion occupation bitmask, bit 2s+c */
  ferm: Uint32Array;
  /** base-3 link code -> first state index with that link configuration */
  linkCodeStart: Map<number, number>;
}

export const MAX_BASIS = 4_000_000;

const POW3: number[] = [];
for (let i = 0, p = 1; i <= 33; i++, p *= 3) POW3.push(p);

export function linkCode(E: ArrayLike<number>, offset: number, nLinks: number): number {
  let code = 0;
  for (let l = 0; l < nLinks; l++) code += (E[offset + l] + 1) * POW3[l];
  return code;
}

export function buildBasis(lattice: Lattice, staticCharges?: ArrayLike<number>): Basis {
  const { nSites, nLinks, linkFrom, linkTo } = lattice;
  if (nLinks > 33) throw new Error('too many links for base-3 link codes');
  const Q = new Int8Array(nSites);
  if (staticCharges) Q.set(staticCharges);

  // a site's charge is known once its highest-index link has been assigned
  const lastLink = new Int32Array(nSites).fill(-1);
  for (let l = 0; l < nLinks; l++) {
    lastLink[linkFrom[l]] = Math.max(lastLink[linkFrom[l]], l);
    lastLink[linkTo[l]] = Math.max(lastLink[linkTo[l]], l);
  }
  const closesAt: number[][] = Array.from({ length: nLinks }, () => []);
  const isolated: number[] = [];
  for (let s = 0; s < nSites; s++) {
    if (lastLink[s] >= 0) closesAt[lastLink[s]].push(s);
    else isolated.push(s);
  }
  for (const s of isolated) if (Q[s] < -1 || Q[s] > 1) return emptyBasis(lattice, Q);

  const linksOut: number[] = [];
  const fermOut: number[] = [];
  const linkCodeStart = new Map<number, number>();
  const E = new Int8Array(nLinks);
  const div = new Int32Array(nSites);

  const emit = () => {
    let zeros = 0;
    const zeroSites: number[] = [];
    let baseFerm = 0;
    for (let s = 0; s < nSites; s++) {
      const qs = lastLink[s] >= 0 ? div[s] - Q[s] : -Q[s];
      if (qs === 1) baseFerm |= 3 << (2 * s);
      else if (qs === 0) {
        zeroSites.push(s);
        zeros++;
      }
    }
    const count = 1 << zeros;
    if (fermOut.length + count > MAX_BASIS) throw new Error('basis too large');
    linkCodeStart.set(linkCode(E, 0, nLinks), fermOut.length);
    for (let r = 0; r < count; r++) {
      let f = baseFerm;
      for (let k = 0; k < zeros; k++) {
        const s = zeroSites[k];
        // bit k of the rank says the upper component (mode 2s) is occupied
        f |= ((r >> k) & 1 ? 1 : 2) << (2 * s);
      }
      fermOut.push(f >>> 0);
      for (let l = 0; l < nLinks; l++) linksOut.push(E[l]);
    }
  };

  const dfs = (l: number) => {
    if (l === nLinks) {
      emit();
      return;
    }
    const a = linkFrom[l];
    const b = linkTo[l];
    for (let e = -1; e <= 1; e++) {
      E[l] = e;
      div[a] += e;
      div[b] -= e;
      let ok = true;
      for (const s of closesAt[l]) {
        const qs = div[s] - Q[s];
        if (qs < -1 || qs > 1) {
          ok = false;
          break;
        }
      }
      if (ok) dfs(l + 1);
      div[a] -= e;
      div[b] += e;
    }
    E[l] = 0;
  };
  dfs(0);

  return {
    lattice,
    staticCharges: Q,
    dim: fermOut.length,
    links: Int8Array.from(linksOut),
    ferm: Uint32Array.from(fermOut),
    linkCodeStart,
  };
}

function emptyBasis(lattice: Lattice, Q: Int8Array): Basis {
  return {
    lattice,
    staticCharges: Q,
    dim: 0,
    links: new Int8Array(0),
    ferm: new Uint32Array(0),
    linkCodeStart: new Map(),
  };
}

/** Index of the state with link values E (length nLinks) and fermion mask f, or -1. */
export function findState(basis: Basis, E: ArrayLike<number>, eOffset: number, f: number): number {
  const { nLinks, nSites } = basis.lattice;
  const start = basis.linkCodeStart.get(linkCode(E, eOffset, nLinks));
  if (start === undefined) return -1;
  let rank = 0;
  let k = 0;
  for (let s = 0; s < nSites; s++) {
    const pair = (f >>> (2 * s)) & 3;
    if (pair === 1) rank |= 1 << k;
    if (pair === 1 || pair === 2) k++;
  }
  const idx = start + rank;
  if (idx >= basis.dim || basis.ferm[idx] !== f >>> 0) return -1;
  return linkCode(basis.links, idx * nLinks, nLinks) === linkCode(E, eOffset, nLinks) ? idx : -1;
}

export const siteCharge = (f: number, s: number): number =>
  (((f >>> (2 * s)) & 1) + ((f >>> (2 * s + 1)) & 1)) - 1;
