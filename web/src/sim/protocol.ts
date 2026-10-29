import type { Params } from '../physics/hamiltonian';
import type { LatticeSpec } from '../physics/lattice';
import type { Observables } from '../physics/observables';

export type Initial =
  /** lowest-energy state of the current sector */
  | { kind: 'ground' }
  /** empty space: no particles, no field */
  | { kind: 'bare' }
  /** flux string between the static charges at `from` (+) and `to` (-) */
  | { kind: 'string'; from: number; to: number };

export interface RunRequest {
  type: 'run';
  id: number;
  lattice: LatticeSpec;
  charges: number[];
  params: Params;
  initial: Initial;
  /** evolve if > 0, otherwise just report the initial state */
  tMax: number;
  dt: number;
}

export interface SampleRequest {
  type: 'sample';
  id: number;
  frame: number;
  count: number;
}

export type Request = RunRequest | SampleRequest;

export interface FrameMsg {
  type: 'frame';
  id: number;
  index: number;
  t: number;
  obs: Observables;
}

export interface ReadyMsg {
  type: 'ready';
  id: number;
  dim: number;
  buildMs: number;
  frames: number;
}

export interface Snapshot {
  E: Int8Array;
  ferm: number;
  probability: number;
}

export interface SampleMsg {
  type: 'sample';
  id: number;
  snapshots: Snapshot[];
}

export interface DoneMsg {
  type: 'done';
  id: number;
  ms: number;
}

export interface ErrorMsg {
  type: 'error';
  id: number;
  message: string;
}

export type Response = FrameMsg | ReadyMsg | SampleMsg | DoneMsg | ErrorMsg;
