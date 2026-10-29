import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { EffectComposer } from 'three/examples/jsm/postprocessing/EffectComposer.js';
import { OutputPass } from 'three/examples/jsm/postprocessing/OutputPass.js';
import { RenderPass } from 'three/examples/jsm/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/examples/jsm/postprocessing/UnrealBloomPass.js';
import type { Lattice } from '../physics/lattice';
import type { Observables } from '../physics/observables';
import type { Snapshot } from '../sim/protocol';
import { COLORS } from './colors';

export type ViewMode = 'flat' | 'torus';

const SPACING = 2.2;
const SPARKS_PER_LINK = 10;

/** A path through lattice coordinates (u, v), possibly broken where it wraps. */
type Segment = [number, number, number, number];

function glowTexture(): THREE.Texture {
  const c = document.createElement('canvas');
  c.width = c.height = 128;
  const g = c.getContext('2d')!;
  const grad = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  grad.addColorStop(0, 'rgba(255,255,255,1)');
  grad.addColorStop(0.2, 'rgba(255,255,255,0.75)');
  grad.addColorStop(0.5, 'rgba(255,255,255,0.18)');
  grad.addColorStop(1, 'rgba(255,255,255,0)');
  g.fillStyle = grad;
  g.fillRect(0, 0, 128, 128);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

interface LinkVis {
  segments: Segment[];
  tubes: THREE.Mesh<THREE.TubeGeometry, THREE.MeshBasicMaterial>[];
  phase: Float64Array;
}

export class LatticeScene {
  readonly renderer: THREE.WebGLRenderer;
  private scene = new THREE.Scene();
  private camera: THREE.PerspectiveCamera;
  private controls: OrbitControls;
  private composer: EffectComposer;
  private bloom: UnrealBloomPass;
  private root = new THREE.Group();
  private snapshotGroup = new THREE.Group();
  private glow = glowTexture();
  private clock = new THREE.Clock();

  private lattice: Lattice | null = null;
  private charges = new Int8Array(0);
  private view: ViewMode = 'flat';
  private obs: Observables | null = null;
  private snapshotFade = 0;

  private siteCores: THREE.Mesh[] = [];
  private pickTargets: THREE.Mesh[] = [];
  private particleSprites: THREE.Sprite[] = [];
  private antiSprites: THREE.Sprite[] = [];
  private staticMeshes: THREE.Object3D[] = [];
  private links: LinkVis[] = [];
  private plaquettes: THREE.Mesh<THREE.BufferGeometry, THREE.MeshBasicMaterial>[] = [];
  private sparks!: THREE.Points;
  private sparkPos!: Float32Array;
  private sparkCol!: Float32Array;

  private raycaster = new THREE.Raycaster();
  private hovered = -1;
  onSiteClick: ((site: number) => void) | null = null;
  onSiteHover: ((site: number, screen: { x: number; y: number } | null) => void) | null = null;

  constructor(private canvas: HTMLCanvasElement) {
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, powerPreference: 'high-performance' });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.scene.background = new THREE.Color(COLORS.space);
    this.scene.fog = new THREE.FogExp2(COLORS.space, 0.012);

    this.camera = new THREE.PerspectiveCamera(45, 1, 0.1, 500);
    this.camera.position.set(0, -9, 9);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.autoRotate = true;
    this.controls.autoRotateSpeed = 0.35;
    this.controls.minDistance = 3;
    this.controls.maxDistance = 60;

    this.composer = new EffectComposer(this.renderer);
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.bloom = new UnrealBloomPass(new THREE.Vector2(1, 1), 0.75, 0.5, 0.22);
    this.composer.addPass(this.bloom);
    this.composer.addPass(new OutputPass());

    this.scene.add(this.root, this.snapshotGroup, this.starfield());
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.6));
    const key = new THREE.PointLight(0xffffff, 60, 0, 1.6);
    key.position.set(4, -6, 10);
    this.scene.add(key);

    canvas.addEventListener('pointermove', (e) => this.pointer(e, false));
    canvas.addEventListener('click', (e) => this.pointer(e, true));
    canvas.addEventListener('pointerdown', () => (this.controls.autoRotate = false));

    new ResizeObserver(() => this.resize()).observe(canvas);
    this.resize();
    this.renderer.setAnimationLoop(() => this.frame());
  }

  setAutoRotate(on: boolean) {
    this.controls.autoRotate = on;
  }

  private starfield(): THREE.Points {
    const n = 2500;
    const pos = new Float32Array(3 * n);
    const rnd = (a: number) => Math.sin(a * 12.9898) * 43758.5453 % 1;
    for (let i = 0; i < n; i++) {
      const u = Math.abs(rnd(i + 1));
      const v = Math.abs(rnd(i + 7.31));
      const th = 2 * Math.PI * u;
      const ph = Math.acos(2 * v - 1);
      const r = 80 + 60 * Math.abs(rnd(i + 3.7));
      pos.set([r * Math.sin(ph) * Math.cos(th), r * Math.sin(ph) * Math.sin(th), r * Math.cos(ph)], 3 * i);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    const m = new THREE.PointsMaterial({ color: 0x8890b0, size: 0.35, sizeAttenuation: true, fog: false });
    return new THREE.Points(g, m);
  }

  private resize() {
    const w = this.canvas.clientWidth;
    const h = this.canvas.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.composer.setSize(w, h);
    this.bloom.setSize(w, h);
    this.camera.aspect = w / h;
    // keep the lattice in the space left free by the side panel and dock
    if (w > 860) this.camera.setViewOffset(w, h, -Math.min(200, w * 0.11), Math.min(60, h * 0.05), w, h);
    else this.camera.setViewOffset(w, h, 0, h * 0.06, w, h);
    this.camera.updateProjectionMatrix();
  }

  // --- geometry -----------------------------------------------------------

  private torusRadii(): [number, number] {
    const L = this.lattice!;
    const R = Math.max(2.4, (L.Lx * SPACING) / (2 * Math.PI) + 0.9);
    const r = L.Ly > 1 ? Math.max(1.1, (L.Ly * SPACING) / (2 * Math.PI)) : 0;
    return [R, Math.min(r, R * 0.7)];
  }

  /** Lattice coordinates -> world position. */
  private map(u: number, v: number, lift = 0): THREE.Vector3 {
    const L = this.lattice!;
    if (this.view === 'torus' && L.periodic) {
      const [R, r] = this.torusRadii();
      const th = (2 * Math.PI * u) / L.Lx;
      const ph = L.Ly > 1 ? (2 * Math.PI * v) / L.Ly + Math.PI / 2 : 0;
      const rr = r + lift;
      return new THREE.Vector3((R + rr * Math.cos(ph)) * Math.cos(th), (R + rr * Math.cos(ph)) * Math.sin(th), rr * Math.sin(ph));
    }
    return new THREE.Vector3((u - (L.Lx - 1) / 2) * SPACING, (v - (L.Ly - 1) / 2) * SPACING, lift);
  }

  private linkSegments(l: number): Segment[] {
    const L = this.lattice!;
    const s = L.linkFrom[l];
    const t = L.linkTo[l];
    const d = L.linkDir[l];
    const [u, v] = [s % L.Lx, Math.floor(s / L.Lx)];
    const [tu, tv] = [t % L.Lx, Math.floor(t / L.Lx)];
    const du = d === 0 ? 1 : 0;
    const dv = d === 1 ? 1 : 0;
    const wraps = tu !== u + du || tv !== v + dv;
    if (!wraps || this.view === 'torus') return [[u, v, u + du, v + dv]];
    return [
      [u, v, u + 0.42 * du, v + 0.42 * dv],
      [tu - 0.42 * du, tv - 0.42 * dv, tu, tv],
    ];
  }

  /** Point at fraction f in [0, 1] along a link's (possibly broken) path. */
  private linkPoint(segs: Segment[], f: number): THREE.Vector3 {
    const seg = segs.length === 1 ? segs[0] : f < 0.5 ? segs[0] : segs[1];
    const g = segs.length === 1 ? f : f < 0.5 ? f * 2 : (f - 0.5) * 2;
    return this.map(seg[0] + (seg[2] - seg[0]) * g, seg[1] + (seg[3] - seg[1]) * g);
  }

  setLattice(L: Lattice, charges: Int8Array, view: ViewMode) {
    this.lattice = L;
    this.charges = charges;
    this.view = L.periodic ? view : 'flat';
    this.obs = null;
    this.rebuild();
  }

  setView(view: ViewMode) {
    if (!this.lattice) return;
    this.view = this.lattice.periodic ? view : 'flat';
    this.rebuild();
  }

  private clear(group: THREE.Group) {
    group.traverse((o) => {
      const m = o as THREE.Mesh;
      if (m.geometry) m.geometry.dispose();
      const mat = m.material as THREE.Material | THREE.Material[] | undefined;
      if (Array.isArray(mat)) mat.forEach((x) => x.dispose());
      else mat?.dispose();
    });
    group.clear();
  }

  private rebuild() {
    const L = this.lattice!;
    this.clear(this.root);
    this.clear(this.snapshotGroup);
    this.siteCores = [];
    this.pickTargets = [];
    this.particleSprites = [];
    this.antiSprites = [];
    this.staticMeshes = [];
    this.links = [];
    this.plaquettes = [];

    const coreGeo = new THREE.SphereGeometry(0.13, 24, 16);
    const pickGeo = new THREE.SphereGeometry(0.55, 8, 6);
    for (let s = 0; s < L.nSites; s++) {
      const p = this.map(s % L.Lx, Math.floor(s / L.Lx));
      const core = new THREE.Mesh(coreGeo, new THREE.MeshBasicMaterial({ color: COLORS.siteIdle }));
      core.position.copy(p);
      this.root.add(core);
      this.siteCores.push(core);

      const pick = new THREE.Mesh(pickGeo, new THREE.MeshBasicMaterial({ visible: false }));
      pick.position.copy(p);
      pick.userData.site = s;
      this.root.add(pick);
      this.pickTargets.push(pick);

      for (const [color, list] of [
        [COLORS.particle, this.particleSprites],
        [COLORS.antiparticle, this.antiSprites],
      ] as const) {
        const sp = new THREE.Sprite(
          new THREE.SpriteMaterial({ map: this.glow, color, blending: THREE.AdditiveBlending, depthWrite: false, transparent: true }),
        );
        sp.position.copy(p);
        sp.scale.setScalar(0.001);
        this.root.add(sp);
        list.push(sp);
      }

      const q = this.charges[s];
      if (q !== 0) this.addStaticCharge(p, q);
    }

    for (let l = 0; l < L.nLinks; l++) {
      const segments = this.linkSegments(l);
      const tubes = segments.map((seg) => {
        const pts = Array.from({ length: 16 }, (_, i) => this.map(seg[0] + ((seg[2] - seg[0]) * i) / 15, seg[1] + ((seg[3] - seg[1]) * i) / 15));
        const curve = new THREE.CatmullRomCurve3(pts);
        const tube = new THREE.Mesh(
          new THREE.TubeGeometry(curve, 24, 0.045, 8, false),
          new THREE.MeshBasicMaterial({ color: COLORS.linkIdle, transparent: true, opacity: 0.9 }),
        );
        this.root.add(tube);
        return tube;
      });
      const phase = new Float64Array(SPARKS_PER_LINK).map((_, i) => i / SPARKS_PER_LINK);
      this.links.push({ segments, tubes, phase });
    }

    for (let p = 0; p < L.nPlaquettes; p++) {
      const s = L.plaquetteSite[p];
      const [u, v] = [s % L.Lx, Math.floor(s / L.Lx)];
      const n = 10;
      const inset = this.view === 'torus' ? 0.3 : 0.12;
      const pos: number[] = [];
      const idx: number[] = [];
      for (let i = 0; i <= n; i++)
        for (let j = 0; j <= n; j++) {
          const q = this.map(u + inset + ((1 - 2 * inset) * i) / n, v + inset + ((1 - 2 * inset) * j) / n, -0.02);
          pos.push(q.x, q.y, q.z);
        }
      for (let i = 0; i < n; i++)
        for (let j = 0; j < n; j++) {
          const a = i * (n + 1) + j;
          idx.push(a, a + n + 1, a + 1, a + 1, a + n + 1, a + n + 2);
        }
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
      g.setIndex(idx);
      const mesh = new THREE.Mesh(
        g,
        new THREE.MeshBasicMaterial({
          color: COLORS.magnetic,
          transparent: true,
          opacity: 0,
          blending: THREE.AdditiveBlending,
          depthWrite: false,
          side: THREE.DoubleSide,
        }),
      );
      this.root.add(mesh);
      this.plaquettes.push(mesh);
    }

    const nSparks = L.nLinks * SPARKS_PER_LINK;
    this.sparkPos = new Float32Array(3 * nSparks);
    this.sparkCol = new Float32Array(3 * nSparks);
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(this.sparkPos, 3));
    g.setAttribute('color', new THREE.BufferAttribute(this.sparkCol, 3));
    this.sparks = new THREE.Points(
      g,
      new THREE.PointsMaterial({
        size: 0.32,
        map: this.glow,
        vertexColors: true,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        transparent: true,
      }),
    );
    this.root.add(this.sparks);
    this.fitCamera();
  }

  private addStaticCharge(p: THREE.Vector3, q: number) {
    const color = q > 0 ? COLORS.particle : COLORS.antiparticle;
    const group = new THREE.Group();
    group.position.copy(p);
    const gem = new THREE.Mesh(
      new THREE.OctahedronGeometry(0.34),
      new THREE.MeshStandardMaterial({ color, emissive: color, emissiveIntensity: 0.35, metalness: 0.4, roughness: 0.3 }),
    );
    const ring = new THREE.Mesh(
      new THREE.TorusGeometry(0.55, 0.025, 8, 48),
      new THREE.MeshBasicMaterial({ color, transparent: true, opacity: 0.8 }),
    );
    group.add(gem, ring);
    this.root.add(group);
    this.staticMeshes.push(group);
  }

  private fitCamera() {
    const box = new THREE.Box3().setFromObject(this.root);
    const radius = box.getSize(new THREE.Vector3()).length() / 2 + 0.8;
    const center = box.getCenter(new THREE.Vector3());
    this.controls.target.copy(center);
    const fov = (this.camera.fov * Math.PI) / 180;
    const fit = radius / Math.sin(fov / 2);
    const narrow = Math.max(1, 1.3 / this.camera.aspect);
    const dist = Math.max(6, fit * 0.95 * narrow);
    const dir = this.view === 'torus' ? new THREE.Vector3(0, -0.75, 0.66) : new THREE.Vector3(0, -0.62, 0.78);
    this.camera.position.copy(center).addScaledVector(dir.normalize(), dist);
    this.controls.update();
  }

  // --- data ---------------------------------------------------------------

  setObservables(obs: Observables | null) {
    this.obs = obs;
  }

  /** Show a single measured configuration on top of the expectation values. */
  showSnapshot(snap: Snapshot | null) {
    this.clear(this.snapshotGroup);
    if (!snap || !this.lattice) {
      this.snapshotFade = 0;
      return;
    }
    const L = this.lattice;
    this.snapshotFade = 1;
    const ball = new THREE.SphereGeometry(0.22, 24, 16);
    for (let s = 0; s < L.nSites; s++) {
      const up = (snap.ferm >>> (2 * s)) & 1;
      const down = (snap.ferm >>> (2 * s + 1)) & 1;
      const p = this.map(s % L.Lx, Math.floor(s / L.Lx));
      const put = (color: THREE.ColorRepresentation, offset: number) => {
        const m = new THREE.Mesh(ball, new THREE.MeshBasicMaterial({ color, transparent: true }));
        m.position.copy(p).add(new THREE.Vector3(offset, 0, 0.35));
        this.snapshotGroup.add(m);
      };
      const hasParticle = up === 1;
      const hasAnti = down === 0;
      if (hasParticle && hasAnti) {
        put(COLORS.particleHot, -0.22);
        put(COLORS.antiparticleHot, 0.22);
      } else if (hasParticle) put(COLORS.particleHot, 0);
      else if (hasAnti) put(COLORS.antiparticleHot, 0);
    }
    const cone = new THREE.ConeGeometry(0.16, 0.42, 16);
    for (let l = 0; l < L.nLinks; l++) {
      const e = snap.E[l];
      if (e === 0) continue;
      const segs = this.linkSegments(l);
      const mid = this.linkPoint(segs, 0.5);
      const ahead = this.linkPoint(segs, e > 0 ? 0.56 : 0.44);
      const m = new THREE.Mesh(cone, new THREE.MeshBasicMaterial({ color: COLORS.fieldHot, transparent: true }));
      m.position.copy(mid);
      m.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), ahead.sub(mid).normalize());
      this.snapshotGroup.add(m);
    }
  }

  // --- interaction ----------------------------------------------------------

  private pointer(e: PointerEvent | MouseEvent, click: boolean) {
    const rect = this.canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const hit = this.raycaster.intersectObjects(this.pickTargets)[0];
    const site = hit ? (hit.object.userData.site as number) : -1;
    if (click && site >= 0) this.onSiteClick?.(site);
    if (site !== this.hovered || !click) {
      this.hovered = site;
      this.canvas.style.cursor = site >= 0 && this.onSiteClick ? 'pointer' : '';
      this.onSiteHover?.(site, site >= 0 ? { x: e.clientX, y: e.clientY } : null);
    }
  }

  // --- animation ----------------------------------------------------------

  private frame() {
    const dt = Math.min(this.clock.getDelta(), 0.1);
    const time = this.clock.elapsedTime;
    this.controls.update();
    if (this.lattice) this.animate(dt, time);
    this.composer.render();
  }

  private animate(dt: number, time: number) {
    const L = this.lattice!;
    const obs = this.obs;
    const fade = 1 - 0.75 * this.snapshotFade;

    for (let s = 0; s < L.nSites; s++) {
      const np = obs ? obs.particles[s] : 0;
      const na = obs ? obs.antiparticles[s] : 0;
      const wobble = 1 + 0.06 * Math.sin(time * 3 + s * 1.7);
      this.particleSprites[s].scale.setScalar((0.3 + 2.6 * Math.sqrt(np)) * wobble * (np > 0.002 ? 1 : 0.001));
      this.antiSprites[s].scale.setScalar((0.3 + 2.6 * Math.sqrt(na)) * wobble * (na > 0.002 ? 1 : 0.001));
      (this.particleSprites[s].material as THREE.SpriteMaterial).opacity = Math.min(1, 0.3 + 1.5 * np) * fade;
      (this.antiSprites[s].material as THREE.SpriteMaterial).opacity = Math.min(1, 0.3 + 1.5 * na) * fade;
      const core = this.siteCores[s].material as THREE.MeshBasicMaterial;
      core.color.set(s === this.hovered ? COLORS.siteHover : COLORS.siteIdle);
    }
    for (const m of this.staticMeshes) {
      m.children[0].rotation.z += dt * 0.8;
      m.children[0].rotation.x += dt * 0.5;
      m.children[1].rotation.x = Math.sin(time * 0.9) * 0.6;
      m.children[1].rotation.y += dt * 1.2;
    }

    const idle = new THREE.Color(COLORS.linkIdle);
    const hot = new THREE.Color(COLORS.field);
    const sparkColor = new THREE.Color(COLORS.fieldHot);
    for (let l = 0; l < L.nLinks; l++) {
      const E = obs ? obs.E[l] : 0;
      const E2 = obs ? obs.E2[l] : 0;
      const vis = this.links[l];
      const c = idle.clone().lerp(hot, Math.min(1, E2));
      for (const t of vis.tubes) t.material.color.copy(c).multiplyScalar(0.35 + 0.75 * E2 * fade);
      const bright = Math.min(1, 0.12 * Math.sqrt(E2) + 0.8 * Math.abs(E)) * fade;
      for (let k = 0; k < SPARKS_PER_LINK; k++) {
        // sparks flow along the field, jitter where the field is only fluctuating
        vis.phase[k] = (((vis.phase[k] + dt * 0.9 * E) % 1) + 1) % 1;
        const jitter = 0.03 * Math.sqrt(Math.max(0, E2 - E * E)) * Math.sin(time * 7 + k * 2.1 + l);
        const f = Math.min(0.97, Math.max(0.03, vis.phase[k] + jitter));
        const p = this.linkPoint(vis.segments, f);
        const i = l * SPARKS_PER_LINK + k;
        this.sparkPos.set([p.x, p.y, p.z], 3 * i);
        const tail = 0.55 + 0.45 * Math.sin(Math.PI * f);
        this.sparkCol.set([sparkColor.r * bright * tail, sparkColor.g * bright * tail, sparkColor.b * bright * tail], 3 * i);
      }
    }
    this.sparks.geometry.attributes.position.needsUpdate = true;
    this.sparks.geometry.attributes.color.needsUpdate = true;

    for (let p = 0; p < L.nPlaquettes; p++) {
      const u = obs ? obs.plaquette[p] : 0;
      const cap = this.view === 'torus' ? 0.25 : 0.4;
      this.plaquettes[p].material.opacity = Math.min(cap, 0.6 * Math.abs(u)) * (0.8 + 0.2 * Math.sin(time * 2 + p)) * fade;
    }

    if (this.snapshotFade > 0) {
      this.snapshotFade = Math.max(0, this.snapshotFade - dt / 2.2);
      const a = Math.min(1, this.snapshotFade * 1.6);
      this.snapshotGroup.traverse((o) => {
        const m = (o as THREE.Mesh).material as THREE.MeshBasicMaterial | undefined;
        if (m) m.opacity = a;
      });
      if (this.snapshotFade === 0) this.clear(this.snapshotGroup);
    }
  }
}
