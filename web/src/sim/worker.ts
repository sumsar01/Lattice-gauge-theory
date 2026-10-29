/// <reference lib="webworker" />
import { buildBasis, findState } from '../physics/basis';
import { applyH, buildHamiltonian, type HamiltonianParts, type Params } from '../physics/hamiltonian';
import { buildLattice } from '../physics/lattice';
import { type CVec, groundState, krylovStep, mulberry32 } from '../physics/linalg';
import { basisVector, diracSea, measure, sample, stringState } from '../physics/observables';
import type { Request, Response, RunRequest } from './protocol';

const ctx = self as unknown as DedicatedWorkerGlobalScope;
const post = (m: Response) => ctx.postMessage(m);

const STATE_BUDGET_BYTES = 256 * 1024 * 1024;

let cachedKey = '';
let cachedH: HamiltonianParts | null = null;
let currentId = 0;
let frames: CVec[] = [];
let framesStride = 1;
const rand = mulberry32(Date.now());

function hamiltonianFor(req: RunRequest): HamiltonianParts {
  const key = JSON.stringify([req.lattice, req.charges]);
  if (key !== cachedKey || !cachedH) {
    cachedH = null;
    const L = buildLattice(req.lattice);
    cachedH = buildHamiltonian(buildBasis(L, req.charges));
    cachedKey = key;
  }
  return cachedH;
}

let trial = new Float64Array(0);

function evolve(apply: (x: CVec, y: CVec) => void, x: CVec, dt: number, depth = 0) {
  if (trial.length !== x.length) trial = new Float64Array(x.length);
  trial.set(x);
  const err = krylovStep(apply, trial, dt, 20);
  if (err < 1e-9 || depth > 6) {
    x.set(trial);
    return;
  }
  evolve(apply, x, dt / 2, depth + 1);
  evolve(apply, x, dt / 2, depth + 1);
}

const yieldToEvents = () => new Promise<void>((r) => setTimeout(r, 0));

async function run(req: RunRequest) {
  // let newer requests queued behind this one supersede it before the expensive part
  await yieldToEvents();
  if (currentId !== req.id) return;
  const t0 = performance.now();
  frames = [];
  const H = hamiltonianFor(req);
  if (H.dim === 0) throw new Error('No allowed states: these static charges cannot be screened.');
  const p: Params = req.params;
  const apply = (x: CVec, y: CVec) => applyH(H, p, x, y);

  let x: CVec;
  if (req.initial.kind === 'ground') {
    x = groundState(apply, H.dim).state;
    if (currentId !== req.id) return;
  }
  else if (req.initial.kind === 'string') x = stringState(H.basis, req.initial.from, req.initial.to);
  else {
    const L = H.basis.lattice;
    const idx = findState(H.basis, new Int8Array(L.nLinks), 0, diracSea(L.nSites));
    if (idx < 0) throw new Error('Bare vacuum is not allowed with static charges present.');
    x = basisVector(H.dim, idx);
  }

  const nFrames = req.tMax > 0 ? Math.round(req.tMax / req.dt) + 1 : 1;
  const keepStates = Math.floor(STATE_BUDGET_BYTES / (16 * H.dim));
  const stride = Math.max(1, Math.ceil(nFrames / Math.max(1, keepStates)));
  framesStride = stride;
  post({ type: 'ready', id: req.id, dim: H.dim, buildMs: performance.now() - t0, frames: nFrames });

  let lastYield = performance.now();
  for (let k = 0; k < nFrames; k++) {
    if (currentId !== req.id) return;
    if (k > 0) evolve(apply, x, req.dt);
    frames[Math.floor(k / stride)] ??= x.slice();
    post({ type: 'frame', id: req.id, index: k, t: k * req.dt, obs: measure(H, p, x) });
    if (performance.now() - lastYield > 25) {
      await yieldToEvents();
      lastYield = performance.now();
      if (currentId !== req.id) return;
    }
  }
  post({ type: 'done', id: req.id, ms: performance.now() - t0 });
}


ctx.onmessage = async (ev: MessageEvent<Request>) => {
  const req = ev.data;
  if (req.type === 'run') {
    currentId = req.id;
    try {
      await run(req);
    } catch (e) {
      post({ type: 'error', id: req.id, message: e instanceof Error ? e.message : String(e) });
    }
  } else if (req.type === 'sample' && cachedH) {
    const k = Math.min(frames.length - 1, Math.round(req.frame / framesStride));
    const x = frames[Math.max(0, k)];
    if (!x) return;
    const snapshots = Array.from({ length: req.count }, () => sample(cachedH!.basis, x, rand()));
    post({ type: 'sample', id: req.id, snapshots });
  }
};
