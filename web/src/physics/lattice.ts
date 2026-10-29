export interface LatticeSpec {
  Lx: number;
  Ly: number;
  periodic: boolean;
}

/** Geometry of a 2D square lattice: sites, directed links (+x, +y) and plaquettes. */
export interface Lattice extends LatticeSpec {
  nSites: number;
  nLinks: number;
  linkFrom: Int32Array;
  linkTo: Int32Array;
  linkDir: Int8Array;
  /** link index leaving site s in direction d at [2s+d], or -1 */
  outLink: Int32Array;
  /** link index arriving at site s from direction d at [2s+d], or -1 */
  inLink: Int32Array;
  /** per plaquette: [bottom x-link, right y-link, top x-link, left y-link] */
  plaquettes: Int32Array;
  nPlaquettes: number;
  /** lower-left corner site of each plaquette */
  plaquetteSite: Int32Array;
}

export const MAX_MODES = 32;

export function siteIndex(L: LatticeSpec, x: number, y: number): number {
  return y * L.Lx + x;
}

export function siteCoords(L: LatticeSpec, s: number): [number, number] {
  return [s % L.Lx, Math.floor(s / L.Lx)];
}

export function buildLattice(spec: LatticeSpec): Lattice {
  const { Lx, Ly, periodic } = spec;
  if (Lx < 2 || Ly < 1) throw new Error('lattice must be at least 2x1');
  if (2 * Lx * Ly > MAX_MODES) throw new Error(`at most ${MAX_MODES / 2} sites supported`);
  const nSites = Lx * Ly;
  const from: number[] = [];
  const to: number[] = [];
  const dir: number[] = [];
  const outLink = new Int32Array(2 * nSites).fill(-1);
  const inLink = new Int32Array(2 * nSites).fill(-1);

  for (let s = 0; s < nSites; s++) {
    const [x, y] = siteCoords(spec, s);
    for (let d = 0; d < 2; d++) {
      let nx = x + (d === 0 ? 1 : 0);
      let ny = y + (d === 1 ? 1 : 0);
      if (d === 1 && Ly === 1) continue;
      if (nx >= Lx || ny >= Ly) {
        if (!periodic) continue;
        nx %= Lx;
        ny %= Ly;
      }
      const t = siteIndex(spec, nx, ny);
      const l = from.length;
      from.push(s);
      to.push(t);
      dir.push(d);
      outLink[2 * s + d] = l;
      inLink[2 * t + d] = l;
    }
  }

  const plaq: number[] = [];
  const plaqSite: number[] = [];
  for (let s = 0; s < nSites; s++) {
    const bottom = outLink[2 * s];
    const left = outLink[2 * s + 1];
    if (bottom < 0 || left < 0) continue;
    const right = outLink[2 * to[bottom] + 1];
    const top = outLink[2 * to[left]];
    if (right < 0 || top < 0) continue;
    plaq.push(bottom, right, top, left);
    plaqSite.push(s);
  }

  return {
    ...spec,
    nSites,
    nLinks: from.length,
    linkFrom: Int32Array.from(from),
    linkTo: Int32Array.from(to),
    linkDir: Int8Array.from(dir),
    outLink,
    inLink,
    plaquettes: Int32Array.from(plaq),
    nPlaquettes: plaqSite.length,
    plaquetteSite: Int32Array.from(plaqSite),
  };
}
