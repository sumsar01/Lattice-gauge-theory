import { findState, type Basis } from './basis';
import { type HamiltonianParts, type Params, massGap, taggedExpectation } from './hamiltonian';
import type { CVec } from './linalg';
import type { Lattice } from './lattice';

export interface Observables {
  /** <q(s)> of the dynamical matter */
  charge: Float64Array;
  /** <number of particles> per site */
  particles: Float64Array;
  /** <number of antiparticles> per site */
  antiparticles: Float64Array;
  /** <E_l> */
  E: Float64Array;
  /** <E_l^2> */
  E2: Float64Array;
  /** <(U_p + U_p^dag)/2> */
  plaquette: Float64Array;
  /** hopping energy per link */
  hopping: Float64Array;
  energy: { mass: number; electric: number; magnetic: number; hopping: number; total: number };
}

export function measure(H: HamiltonianParts, p: Params, x: CVec): Observables {
  const { basis } = H;
  const L = basis.lattice;
  const { nSites, nLinks, nPlaquettes } = L;
  const charge = new Float64Array(nSites);
  const particles = new Float64Array(nSites);
  const antiparticles = new Float64Array(nSites);
  const E = new Float64Array(nLinks);
  const E2 = new Float64Array(nLinks);
  let mass = 0;
  let electric = 0;

  for (let i = 0; i < H.dim; i++) {
    const pr = x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1];
    if (pr < 1e-16) continue;
    const f = basis.ferm[i];
    for (let s = 0; s < nSites; s++) {
      const up = (f >>> (2 * s)) & 1;
      const down = (f >>> (2 * s + 1)) & 1;
      charge[s] += pr * (up + down - 1);
      particles[s] += pr * up;
      antiparticles[s] += pr * (1 - down);
    }
    const off = i * nLinks;
    for (let l = 0; l < nLinks; l++) {
      const e = basis.links[off + l];
      E[l] += pr * e;
      E2[l] += pr * e * e;
    }
    mass += pr * H.diagMass[i];
    electric += pr * H.diagElectric[i];
  }

  const plaquette = taggedExpectation(H.plaq, x, nPlaquettes);
  let magnetic = 0;
  for (let k = 0; k < nPlaquettes; k++) {
    magnetic += plaquette[k];
    plaquette[k] /= 2;
  }
  const hopping = taggedExpectation(H.hopA, x, nLinks);
  if (p.r !== 0) {
    const hr = taggedExpectation(H.hopR, x, nLinks);
    for (let l = 0; l < nLinks; l++) hopping[l] += p.r * hr[l];
  }
  let hop = 0;
  for (let l = 0; l < nLinks; l++) hop += hopping[l];

  const energy = {
    mass: massGap(H, p) * mass,
    electric: ((p.g * p.g) / 2) * electric,
    magnetic: (-1 / (2 * p.g * p.g)) * magnetic,
    hopping: hop,
    total: 0,
  };
  energy.total = energy.mass + energy.electric + energy.magnetic + energy.hopping;
  return { charge, particles, antiparticles, E, E2, plaquette, hopping, energy };
}

/** Fermion mask of the bare Dirac sea: every lower component filled. */
export function diracSea(nSites: number): number {
  let f = 0;
  for (let s = 0; s < nSites; s++) f |= 2 << (2 * s);
  return f >>> 0;
}

/** Normalised basis state as a complex vector. */
export function basisVector(dim: number, idx: number): CVec {
  const v = new Float64Array(2 * dim);
  v[2 * idx] = 1;
  return v;
}

/**
 * Flux string between static charges: walks x first, then y, from the
 * positive to the negative charge and puts one unit of flux on each link.
 */
export function fluxString(L: Lattice, from: number, to: number): Int8Array {
  const E = new Int8Array(L.nLinks);
  let [x, y] = [from % L.Lx, Math.floor(from / L.Lx)];
  const [tx, ty] = [to % L.Lx, Math.floor(to / L.Lx)];
  const step = (target: number, size: number, cur: number): [number, number] => {
    let delta = target - cur;
    if (L.periodic && Math.abs(delta) > size / 2) delta -= Math.sign(delta) * size;
    return [Math.sign(delta), Math.abs(delta)];
  };
  for (const d of [0, 1] as const) {
    const [sgn, n] = d === 0 ? step(tx, L.Lx, x) : step(ty, L.Ly, y);
    for (let k = 0; k < n; k++) {
      const s = y * L.Lx + x;
      if (sgn > 0) {
        E[L.outLink[2 * s + d]] += 1;
        if (d === 0) x = (x + 1) % L.Lx;
        else y = (y + 1) % L.Ly;
      } else {
        E[L.inLink[2 * s + d]] -= 1;
        if (d === 0) x = (x - 1 + L.Lx) % L.Lx;
        else y = (y - 1 + L.Ly) % L.Ly;
      }
    }
  }
  return E;
}

export function stringState(basis: Basis, from: number, to: number): CVec {
  const E = fluxString(basis.lattice, from, to);
  const idx = findState(basis, E, 0, diracSea(basis.lattice.nSites));
  if (idx < 0) throw new Error('string state not in basis');
  return basisVector(basis.dim, idx);
}

/** Draw one classical configuration from |psi|^2, as if we looked. */
export function sample(basis: Basis, x: CVec, u: number): { E: Int8Array; ferm: number; probability: number } {
  let total = 0;
  for (let i = 0; i < 2 * basis.dim; i++) total += x[i] * x[i];
  u *= total;
  let acc = 0;
  let idx = basis.dim - 1;
  for (let i = 0; i < basis.dim; i++) {
    acc += x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1];
    if (acc >= u) {
      idx = i;
      break;
    }
  }
  const n = basis.lattice.nLinks;
  return {
    E: basis.links.slice(idx * n, (idx + 1) * n),
    ferm: basis.ferm[idx],
    probability: x[2 * idx] ** 2 + x[2 * idx + 1] ** 2,
  };
}
