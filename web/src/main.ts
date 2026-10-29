import { buildLattice, type Lattice } from './physics/lattice';
import type { Params } from './physics/hamiltonian';
import type { Observables } from './physics/observables';
import type { Response, RunRequest, Snapshot } from './sim/protocol';
import { SCENARIOS, type Scenario } from './sim/scenarios';
import { TimeChart } from './ui/chart';
import { LatticeScene, type ViewMode } from './viz/scene';

const $ = <T extends HTMLElement = HTMLElement>(id: string) => document.getElementById(id) as T;

let scene: LatticeScene;
try {
  scene = new LatticeScene($('scene') as HTMLCanvasElement);
} catch (e) {
  document.body.innerHTML = `<div class="fatal"><p>This universe needs WebGL, which your browser could not start.<br>${String(e)}</p></div>`;
  throw e;
}
const chart = new TimeChart($('chart'));
const worker = new Worker(new URL('./sim/worker.ts', import.meta.url), { type: 'module' });

interface Frame {
  t: number;
  obs: Observables;
}

const state = {
  scenario: SCENARIOS[0],
  latticeIdx: 0,
  lattice: null as Lattice | null,
  params: { ...SCENARIOS[0].params } as Params,
  charges: [] as number[],
  view: 'torus' as ViewMode,
  runId: 0,
  frames: [] as Frame[],
  totalFrames: 0,
  dim: 0,
  computeMs: 0,
  playing: true,
  playT: 0,
  speed: 1,
};

const sum = (a: Float64Array) => a.reduce((x, y) => x + y, 0);

// --- worker -----------------------------------------------------------------

function run() {
  const sc = state.scenario;
  const spec = sc.lattices[state.latticeIdx].spec;
  state.lattice = buildLattice(spec);
  state.runId++;
  state.frames = [];
  state.totalFrames = 0;
  state.playT = 0;
  scene.setLattice(state.lattice, Int8Array.from(state.charges), state.view);
  scene.showSnapshot(null);
  chart.reset(sc.tMax || 1);
  $('look-caption').textContent = '';
  $('dock').classList.toggle('static', sc.tMax === 0);
  $('view-toggle').style.display = spec.periodic ? '' : 'none';
  status(sc.lattices[state.latticeIdx].heavy ? 'Solving a big universe, this can take a while…' : 'Solving…');
  const req: RunRequest = {
    type: 'run',
    id: state.runId,
    lattice: spec,
    charges: state.charges,
    params: { ...state.params },
    initial: sc.initial(spec, state.charges),
    tMax: sc.tMax,
    dt: sc.dt,
  };
  worker.postMessage(req);
}

worker.onmessage = (ev: MessageEvent<Response>) => {
  const msg = ev.data;
  if (msg.id !== state.runId) return;
  switch (msg.type) {
    case 'ready':
      state.dim = msg.dim;
      state.totalFrames = msg.frames;
      state.computeMs = msg.buildMs;
      status('');
      break;
    case 'frame':
      state.frames[msg.index] = { t: msg.t, obs: msg.obs };
      chart.push(msg.t, sum(msg.obs.particles) + sum(msg.obs.antiparticles), msg.obs.energy.electric);
      if (state.scenario.tMax === 0) {
        scene.setObservables(msg.obs);
        renderTiles(msg.obs);
      }
      break;
    case 'done':
      state.computeMs = msg.ms;
      renderNerdStats();
      break;
    case 'sample':
      showSnapshot(msg.snapshots[0]);
      break;
    case 'error':
      status(msg.message, true);
      break;
  }
};

// --- playback -----------------------------------------------------------------

function lerp(a: Float64Array, b: Float64Array, f: number): Float64Array {
  const out = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) out[i] = a[i] + (b[i] - a[i]) * f;
  return out;
}

function observablesAt(t: number): Observables | null {
  const { frames } = state;
  if (!frames.length) return null;
  const dt = state.scenario.dt;
  const i = Math.min(Math.floor(t / dt), frames.length - 1);
  const a = frames[i];
  const b = frames[i + 1];
  if (!b) return a.obs;
  const f = (t - a.t) / dt;
  const lo = a.obs;
  const hi = b.obs;
  const e = (k: keyof Observables['energy']) => lo.energy[k] + (hi.energy[k] - lo.energy[k]) * f;
  return {
    charge: lerp(lo.charge, hi.charge, f),
    particles: lerp(lo.particles, hi.particles, f),
    antiparticles: lerp(lo.antiparticles, hi.antiparticles, f),
    E: lerp(lo.E, hi.E, f),
    E2: lerp(lo.E2, hi.E2, f),
    plaquette: lerp(lo.plaquette, hi.plaquette, f),
    hopping: lerp(lo.hopping, hi.hopping, f),
    energy: { mass: e('mass'), electric: e('electric'), magnetic: e('magnetic'), hopping: e('hopping'), total: e('total') },
  };
}

let lastTick = performance.now();
let lastTileUpdate = 0;
function tick(now: number) {
  const dtReal = Math.min(0.1, (now - lastTick) / 1000);
  lastTick = now;
  const sc = state.scenario;
  if (sc.tMax > 0 && state.frames.length) {
    const available = state.frames[state.frames.length - 1].t;
    if (state.playing) {
      state.playT = Math.min(state.playT + dtReal * state.speed, available);
      if (state.playT >= sc.tMax - 1e-9) setPlaying(false);
    }
    const buffering = state.playing && state.playT >= available && state.frames.length < state.totalFrames;
    status(buffering ? 'Computing the future…' : '');
    const obs = observablesAt(state.playT);
    scene.setObservables(obs);
    chart.setPlayhead(state.playT);
    ($('scrub') as HTMLInputElement).value = String(state.playT / sc.tMax);
    $('clock').textContent = `t = ${state.playT.toFixed(1)}`;
    if (obs && now - lastTileUpdate > 120) {
      renderTiles(obs);
      lastTileUpdate = now;
    }
  }
  requestAnimationFrame(tick);
}
requestAnimationFrame(tick);

function setPlaying(on: boolean) {
  const sc = state.scenario;
  if (on && state.playT >= sc.tMax - 1e-9) state.playT = 0;
  state.playing = on;
  $('play').textContent = on ? '❚❚' : '▶';
  $('play').setAttribute('aria-label', on ? 'Pause' : 'Play');
}

// --- looking (measurement) -------------------------------------------------------

function look() {
  const frame = Math.round(state.playT / state.scenario.dt);
  worker.postMessage({ type: 'sample', id: state.runId, frame, count: 1 });
}

function showSnapshot(s: Snapshot) {
  scene.showSnapshot(s);
  const L = state.lattice!;
  let p = 0;
  let a = 0;
  for (let k = 0; k < L.nSites; k++) {
    if ((s.ferm >>> (2 * k)) & 1) p++;
    if (!((s.ferm >>> (2 * k + 1)) & 1)) a++;
  }
  const lines = s.E.reduce((n, e) => n + (e !== 0 ? 1 : 0), 0);
  const pct = s.probability * 100;
  const chance = pct >= 0.1 ? `${pct.toFixed(1)}%` : '< 0.1%';
  $('look-caption').innerHTML = `Saw <b>${p}</b> particle${p === 1 ? '' : 's'}, <b>${a}</b> anti, <b>${lines}</b> field line${lines === 1 ? '' : 's'}. Chance: ${chance}`;
}

let keepLooking: number | undefined;
($('keep-looking') as HTMLInputElement).addEventListener('change', (e) => {
  clearInterval(keepLooking);
  if ((e.target as HTMLInputElement).checked) keepLooking = window.setInterval(look, 1400);
});
$('look').addEventListener('click', look);

// --- panels ---------------------------------------------------------------------

function renderTiles(obs: Observables) {
  const n = sum(obs.particles) + sum(obs.antiparticles);
  const fieldShare = obs.energy.electric;
  $('tiles').innerHTML = [
    tile('Particles', n.toFixed(2), '+ antiparticles'),
    tile('Field energy', fieldShare.toFixed(2), 'electric'),
    tile('States', state.dim.toLocaleString(), 'possibilities'),
  ].join('');
}

const tile = (k: string, v: string, s: string) => `<div class="tile"><div class="k">${k}</div><div class="v">${v}</div><div class="s">${s}</div></div>`;

function renderNerdStats() {
  const last = state.frames[0]?.obs;
  const rows: [string, string][] = [
    ['Hilbert space', `${state.dim.toLocaleString()} gauge-invariant states`],
    ['Solve time', `${state.computeMs.toFixed(0)} ms`],
  ];
  if (last) rows.push(['⟨H⟩ at t = 0', last.energy.total.toFixed(5)]);
  $('nerd-stats').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
}

function status(text: string, error = false) {
  const el = $('status');
  if (text) el.textContent = text;
  el.classList.toggle('show', !!text);
  el.classList.toggle('error', error);
}

function selectScenario(sc: Scenario) {
  state.scenario = sc;
  state.latticeIdx = 0;
  state.params = { ...sc.params };
  state.charges = sc.charges(sc.lattices[0].spec);
  document.querySelectorAll<HTMLButtonElement>('#scenarios button').forEach((b) => b.classList.toggle('on', b.dataset.id === sc.id));
  $('title').textContent = sc.title;
  $('tagline').textContent = sc.tagline;
  $('story').innerHTML = sc.story.map((p) => `<p>${p}</p>`).join('');
  const sel = $('lattice') as HTMLSelectElement;
  sel.innerHTML = sc.lattices.map((l, i) => `<option value="${i}">${l.label}</option>`).join('');
  syncSliders();
  scene.onSiteClick = sc.editableCharges ? toggleCharge : null;
  ($('keep-looking') as HTMLInputElement).checked = false;
  clearInterval(keepLooking);
  setPlaying(true);
  run();
}

function toggleCharge(site: number) {
  const q = state.charges[site];
  state.charges[site] = q === 0 ? 1 : q === 1 ? -1 : 0;
  run();
}

function syncSliders() {
  for (const [id, key] of [
    ['g', 'g'],
    ['m', 'mass'],
    ['r', 'r'],
  ] as const) {
    ($(id) as HTMLInputElement).value = String(state.params[key]);
    $(`${id}-val`).textContent = state.params[key].toFixed(2);
  }
}

let debounce: number | undefined;
for (const [id, key] of [
  ['g', 'g'],
  ['m', 'mass'],
  ['r', 'r'],
] as const) {
  $(id).addEventListener('input', (e) => {
    state.params[key] = Number((e.target as HTMLInputElement).value);
    $(`${id}-val`).textContent = state.params[key].toFixed(2);
    clearTimeout(debounce);
    debounce = window.setTimeout(() => {
      setPlaying(true);
      run();
    }, 220);
  });
}

$('lattice').addEventListener('change', (e) => {
  state.latticeIdx = Number((e.target as HTMLSelectElement).value);
  state.charges = state.scenario.charges(state.scenario.lattices[state.latticeIdx].spec);
  setPlaying(true);
  run();
});

$('scenarios').innerHTML = SCENARIOS.map((s) => `<button data-id="${s.id}">${s.title}<small>${s.tagline}</small></button>`).join('');
$('scenarios').addEventListener('click', (e) => {
  const id = (e.target as HTMLElement).closest('button')?.dataset.id;
  const sc = SCENARIOS.find((s) => s.id === id);
  if (sc && sc !== state.scenario) selectScenario(sc);
});

$('play').addEventListener('click', () => setPlaying(!state.playing));
$('scrub').addEventListener('input', (e) => {
  setPlaying(false);
  const avail = state.frames.length ? state.frames[state.frames.length - 1].t : 0;
  state.playT = Math.min(avail, Number((e.target as HTMLInputElement).value) * state.scenario.tMax);
});
$('speed').addEventListener('change', (e) => (state.speed = Number((e.target as HTMLSelectElement).value)));
chart.onSeek = (t) => {
  const avail = state.frames.length ? state.frames[state.frames.length - 1].t : 0;
  state.playT = Math.min(t, avail);
};

$('view-toggle').addEventListener('click', (e) => {
  const v = (e.target as HTMLElement).closest('button')?.dataset.view as ViewMode | undefined;
  if (!v) return;
  state.view = v;
  document.querySelectorAll<HTMLButtonElement>('#view-toggle button').forEach((b) => b.classList.toggle('on', b.dataset.view === v));
  scene.setView(v);
});

const siteTip = $('site-tip');
scene.onSiteHover = (site, at) => {
  const obs = state.scenario.tMax > 0 ? observablesAt(state.playT) : state.frames[0]?.obs;
  if (site < 0 || !at || !obs) {
    siteTip.style.display = 'none';
    return;
  }
  const q = state.charges[site];
  siteTip.innerHTML =
    `particle <b>${(obs.particles[site] * 100).toFixed(0)}%</b> · antiparticle <b>${(obs.antiparticles[site] * 100).toFixed(0)}%</b>` +
    (q ? `<br>fixed charge <b>${q > 0 ? '+1' : '−1'}</b>` : '') +
    (state.scenario.editableCharges ? `<br><span style="color:var(--ink-3)">click to change the fixed charge</span>` : '');
  siteTip.style.display = 'block';
  siteTip.style.left = `${at.x + 14}px`;
  siteTip.style.top = `${at.y + 14}px`;
};

const intro = $('intro');
let seen = false;
try {
  seen = localStorage.getItem('pocket-universe-intro') === '1';
} catch {
  /* storage unavailable */
}
if (seen) intro.classList.add('hidden');
$('intro-go').addEventListener('click', () => {
  intro.classList.add('hidden');
  try {
    localStorage.setItem('pocket-universe-intro', '1');
  } catch {
    /* storage unavailable */
  }
});

selectScenario(SCENARIOS[0]);
