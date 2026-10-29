import { describe, expect, it } from 'vitest';
import { buildBasis, findState } from './basis';
import { applyH, buildHamiltonian, type Params } from './hamiltonian';
import { buildLattice } from './lattice';
import { groundState, krylovStep } from './linalg';
import { basisVector, diracSea } from './observables';

// Values from tools/brute_force_reference.py (full tensor-product space,
// explicit Jordan-Wigner strings, Gauss law imposed by projection).
const CASES: { L: [number, number, boolean]; Q?: number[]; dim: number; E0: number; amp?: [number, number] }[] = [
  { L: [2, 1, false], dim: 6, E0: -2.74816253, amp: [-0.8429682295869357, -0.42860050768865293] },
  { L: [2, 1, true], dim: 16, E0: -2.944342, amp: [-0.6349027438662684, -0.5062104350503419] },
  { L: [3, 1, true], dim: 48, E0: -4.3647471, amp: [0.681736009006989, -0.4495200072068671] },
  { L: [4, 1, false], dim: 68, E0: -5.652977, amp: [0.33983944633783275, 0.7529276925254649] },
  { L: [2, 2, false], dim: 152, E0: -9.65107486, amp: [0.9150905017661157, -0.021238364789779238] },
  { L: [3, 1, false], Q: [1, 0, -1], dim: 14, E0: -2.81613364 },
  { L: [2, 2, true], dim: 6336, E0: -10.34785282, amp: [0.46286527489718055, 0.4933424014511792] },
];
const P: Params = { mass: 0.3, g: 1.2, r: 1 };

describe('matches brute-force reference', () => {
  for (const c of CASES) {
    it(`${c.L[0]}x${c.L[1]} ${c.L[2] ? 'periodic' : 'open'}${c.Q ? ' with static charges' : ''}`, () => {
      const lat = buildLattice({ Lx: c.L[0], Ly: c.L[1], periodic: c.L[2] });
      const B = buildBasis(lat, c.Q);
      expect(B.dim).toBe(c.dim);
      const H = buildHamiltonian(B);
      const apply = (x: Float64Array, y: Float64Array) => applyH(H, P, x, y);
      expect(groundState(apply, H.dim).energy).toBeCloseTo(c.E0, 6);
      if (c.amp) {
        const vac = findState(B, new Int8Array(lat.nLinks), 0, diracSea(lat.nSites));
        const x = basisVector(H.dim, vac);
        for (let k = 0; k < 13; k++) krylovStep(apply, x, 0.1);
        expect(x[2 * vac]).toBeCloseTo(c.amp[0], 7);
        expect(x[2 * vac + 1]).toBeCloseTo(c.amp[1], 7);
      }
    });
  }
});
