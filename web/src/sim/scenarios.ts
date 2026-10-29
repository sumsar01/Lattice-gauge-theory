import type { Params } from '../physics/hamiltonian';
import type { LatticeSpec } from '../physics/lattice';
import type { Initial } from './protocol';

export interface LatticeOption {
  label: string;
  spec: LatticeSpec;
  /** slow lattices get a warning in the UI */
  heavy?: boolean;
}

export interface Scenario {
  id: string;
  title: string;
  tagline: string;
  story: string[];
  lattices: LatticeOption[];
  params: Params;
  /** static background charges for a lattice */
  charges: (L: LatticeSpec) => number[];
  initial: (L: LatticeSpec, charges: number[]) => Initial;
  tMax: number;
  dt: number;
  editableCharges?: boolean;
}

const sites = (L: LatticeSpec) => L.Lx * L.Ly;
const none = (L: LatticeSpec) => new Array(sites(L)).fill(0);

export const SCENARIOS: Scenario[] = [
  {
    id: 'string',
    title: 'Snap the string',
    tagline: 'Why you never find a lonely quark',
    story: [
      'Two heavy charges (the spinning gems) are tied together by a tube of electric field. Its sparks flow from + to −.',
      'Keeping that field alive costs energy. When it gets expensive enough, empty space would rather spend the energy on making a particle–antiparticle pair.',
      'Watch orange and blue clouds bloom out of nothing and cling to the gems, so the string snaps. Quarks are trapped in protons by a stronger version of the same effect.',
    ],
    lattices: [
      { label: '3 × 2 grid', spec: { Lx: 3, Ly: 2, periodic: false } },
      { label: '5 × 1 chain', spec: { Lx: 5, Ly: 1, periodic: false } },
      { label: '4 × 2 grid (slow)', spec: { Lx: 4, Ly: 2, periodic: false }, heavy: true },
    ],
    params: { mass: 0.1, g: 2, r: 0.5 },
    charges: (L) => {
      const q = none(L);
      q[0] = 1;
      q[L.Lx - 1] = -1;
      return q;
    },
    initial: (L) => ({ kind: 'string', from: 0, to: L.Lx - 1 }),
    tMax: 20,
    dt: 0.05,
  },
  {
    id: 'pairs',
    title: 'Something from nothing',
    tagline: 'Empty space boils over',
    story: [
      'We start with perfectly empty space: no particles and no field. Then we let the rules of quantum electrodynamics run.',
      'Pairs of particles (orange) and antiparticles (blue) pop into existence, always together, so the total charge stays zero. Electric field lines (purple sparks) connect each pair.',
      'This universe is a doughnut: walk off the right edge and you come back on the left. Drag to look around it.',
    ],
    lattices: [
      { label: '2 × 2 doughnut', spec: { Lx: 2, Ly: 2, periodic: true } },
      { label: '3 × 2 grid', spec: { Lx: 3, Ly: 2, periodic: false } },
      { label: '8 × 1 ring', spec: { Lx: 8, Ly: 1, periodic: true } },
      { label: '3 × 2 doughnut (slow)', spec: { Lx: 3, Ly: 2, periodic: true }, heavy: true },
    ],
    params: { mass: 0.1, g: 1, r: 0.5 },
    charges: none,
    initial: () => ({ kind: 'bare' }),
    tMax: 20,
    dt: 0.05,
  },
  {
    id: 'vacuum',
    title: 'The restless vacuum',
    tagline: 'Nothing is never quite nothing',
    story: [
      'This is the calmest state the universe can be in, its ground state. Yet the faint clouds show particles and fields are never quite absent.',
      'Press “Look” to measure. Each look gives a different, random snapshot of what is there, and the glow is the average over all possible looks.',
      'Turn the force strength up and the vacuum calms down. Turn it down and the magnetic squares (pink) light up.',
    ],
    lattices: [
      { label: '2 × 2 doughnut', spec: { Lx: 2, Ly: 2, periodic: true } },
      { label: '3 × 2 grid', spec: { Lx: 3, Ly: 2, periodic: false } },
      { label: '8 × 1 ring', spec: { Lx: 8, Ly: 1, periodic: true } },
      { label: '3 × 3 grid (slow)', spec: { Lx: 3, Ly: 3, periodic: false }, heavy: true },
    ],
    params: { mass: 0.1, g: 1, r: 0.5 },
    charges: none,
    initial: () => ({ kind: 'ground' }),
    tMax: 0,
    dt: 0.05,
  },
  {
    id: 'sandbox',
    title: 'Charge sandbox',
    tagline: 'Place charges, watch space react',
    story: [
      'Click a dot to place a heavy + charge, click again for a − charge, and a third time to remove it.',
      'The universe settles into its calmest arrangement around your charges. Sometimes it draws field lines between them. Sometimes it shields them with pairs pulled out of the vacuum.',
      'Gauss’s law always holds: every field line that leaves a charge must end on another charge.',
    ],
    lattices: [
      { label: '3 × 2 grid', spec: { Lx: 3, Ly: 2, periodic: false } },
      { label: '4 × 2 grid (slow)', spec: { Lx: 4, Ly: 2, periodic: false }, heavy: true },
      { label: '2 × 2 doughnut', spec: { Lx: 2, Ly: 2, periodic: true } },
      { label: '6 × 1 chain', spec: { Lx: 6, Ly: 1, periodic: false } },
    ],
    params: { mass: 0.1, g: 1.5, r: 0.5 },
    charges: (L) => {
      const q = none(L);
      q[0] = 1;
      q[sites(L) - 1] = -1;
      return q;
    },
    initial: () => ({ kind: 'ground' }),
    tMax: 0,
    dt: 0.05,
    editableCharges: true,
  },
];
