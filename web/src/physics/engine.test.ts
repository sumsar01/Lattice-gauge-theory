import { describe, expect, it } from 'vitest';
import { buildBasis } from './basis';
import { applyH, buildHamiltonian, type Params } from './hamiltonian';
import { buildLattice } from './lattice';
import { groundState, krylovStep, norm } from './linalg';
import { basisVector, diracSea, measure, stringState } from './observables';
import { findState, siteCharge } from './basis';

const P: Params = { mass: 0.3, g: 1.2, r: 1 };

function setup(Lx: number, Ly: number, periodic: boolean, Q?: number[]) {
  const L = buildLattice({ Lx, Ly, periodic });
  const B = buildBasis(L, Q);
  const H = buildHamiltonian(B);
  const apply = (x: Float64Array, y: Float64Array) => applyH(H, P, x, y);
  return { L, B, H, apply };
}

describe('basis', () => {
  it('counts gauge-invariant states', () => {
    for (const [Lx, Ly, per] of [
      [2, 1, false],
      [2, 2, true],
      [3, 2, false],
      [3, 2, true],
    ] as const) {
      const { B } = setup(Lx, Ly, per);
      expect(B.dim).toBeGreaterThan(0);
      // every state satisfies Gauss law
      const { nSites, nLinks, linkFrom, linkTo } = B.lattice;
      for (let i = 0; i < B.dim; i++) {
        const div = new Int32Array(nSites);
        for (let l = 0; l < nLinks; l++) {
          div[linkFrom[l]] += B.links[i * nLinks + l];
          div[linkTo[l]] -= B.links[i * nLinks + l];
        }
        for (let s = 0; s < nSites; s++) {
          if (div[s] !== siteCharge(B.ferm[i], s)) throw new Error(`Gauss law broken in state ${i}`);
        }
        if (findState(B, B.links, i * nLinks, B.ferm[i]) !== i) throw new Error(`lookup failed for ${i}`);
      }
    }
  });
});

describe('hamiltonian', () => {
  it('is Hermitian', () => {
    const { H, apply } = setup(3, 2, true);
    const n = H.dim;
    const x = new Float64Array(2 * n).map((_, i) => Math.sin(i * 1.3));
    const y = new Float64Array(2 * n).map((_, i) => Math.cos(i * 0.7));
    const Hx = new Float64Array(2 * n);
    const Hy = new Float64Array(2 * n);
    apply(x, Hx);
    apply(y, Hy);
    // <y|Hx> == conj(<x|Hy>)
    let ar = 0, ai = 0, br = 0, bi = 0;
    for (let i = 0; i < 2 * n; i += 2) {
      ar += y[i] * Hx[i] + y[i + 1] * Hx[i + 1];
      ai += y[i] * Hx[i + 1] - y[i + 1] * Hx[i];
      br += x[i] * Hy[i] + x[i + 1] * Hy[i + 1];
      bi += x[i] * Hy[i + 1] - x[i + 1] * Hy[i];
    }
    expect(ar).toBeCloseTo(br, 8);
    expect(ai).toBeCloseTo(-bi, 8);
  });

  it('finds a converged ground state whose measured energy matches', () => {
    const { H, apply } = setup(2, 2, true);
    const gs = groundState(apply, H.dim);
    expect(gs.residual).toBeLessThan(1e-6);
    const obs = measure(H, P, gs.state);
    expect(obs.energy.total).toBeCloseTo(gs.energy, 8);
    for (const q of obs.charge) expect(Math.abs(q)).toBeLessThan(1e-6 + 1);
  });
});

describe('time evolution', () => {
  it('conserves norm and energy', () => {
    const { H, B, apply } = setup(3, 2, false);
    const vac = findState(B, new Int8Array(B.lattice.nLinks), 0, diracSea(B.lattice.nSites));
    const x = basisVector(H.dim, vac);
    const E0 = measure(H, P, x).energy.total;
    for (let k = 0; k < 20; k++) krylovStep(apply, x, 0.1);
    expect(norm(x)).toBeCloseTo(1, 10);
    expect(measure(H, P, x).energy.total).toBeCloseTo(E0, 8);
  });

  it('builds a flux string between static charges', () => {
    const { L, H } = setup(4, 1, false, [1, 0, 0, -1]);
    const x = stringState(H.basis, 0, 3);
    const obs = measure(H, P, x);
    expect(Array.from(obs.E)).toEqual([1, 1, 1]);
    expect(L.nLinks).toBe(3);
  });
});
