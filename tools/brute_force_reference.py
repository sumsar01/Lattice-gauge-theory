"""Brute-force reference: full tensor-product space, explicit JW strings, Gauss projection."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.linalg import expm

def lattice(Lx, Ly, per):
    links = []  # (from, to, dir)
    out = {}
    for s in range(Lx*Ly):
        x, y = s % Lx, s // Lx
        for d in (0, 1):
            if d == 1 and Ly == 1: continue
            nx, ny = x + (d == 0), y + (d == 1)
            if nx >= Lx or ny >= Ly:
                if not per: continue
                nx %= Lx; ny %= Ly
            out[(s, d)] = len(links)
            links.append((s, ny*Lx+nx, d))
    plaq = []
    for s in range(Lx*Ly):
        if (s,0) in out and (s,1) in out:
            b, l = out[(s,0)], out[(s,1)]
            r = out.get((links[b][1],1)); t = out.get((links[l][1],0))
            if r is not None and t is not None: plaq.append((b, r, t, l))
    return links, plaq

def run(Lx, Ly, per, m=0.3, g=1.2, r=1.0, Q=None, t=1.3):
    N = Lx*Ly; links, plaq = lattice(Lx, Ly, per); nm, nl = 2*N, len(links)
    dims = [2]*nm + [3]*nl
    I2, I3 = sp.identity(2, format='csr'), sp.identity(3, format='csr')
    Z = sp.csr_matrix(np.diag([1., -1.]))
    a = sp.csr_matrix(np.array([[0., 1.], [0., 0.]]))  # |0>,|1> ; a|1> = |0>
    Up = sp.csr_matrix(np.diag([1., 1.], -1))  # E=-1,0,1 -> raise
    Esz = sp.csr_matrix(np.diag([-1., 0., 1.]))
    def op(factors):
        mats = [factors.get(k, I2 if k < nm else I3) for k in range(nm+nl)]
        out = mats[0]
        for M in mats[1:]: out = sp.kron(out, M, format='csr')
        return out
    def c(mu):
        f = {k: Z for k in range(mu)}; f[mu] = a; return f
    # Gauss subspace indices
    D = int(np.prod(dims))
    idx = np.arange(D); digits = []
    rem = idx.copy()
    for d in reversed(dims):
        digits.append(rem % d); rem //= d
    digits = digits[::-1]
    n = digits[:nm]; E = [dg - 1 for dg in digits[nm:]]
    Qv = Q or [0]*N
    ok = np.ones(D, bool)
    for s in range(N):
        div = np.zeros(D, int)
        for l, (f_, t_, _) in enumerate(links):
            if f_ == s: div += E[l]
            if t_ == s: div -= E[l]
        ok &= div == n[2*s] + n[2*s+1] - 1 + Qv[s]
    keep = np.nonzero(ok)[0]
    P = sp.csr_matrix((np.ones(len(keep)), (np.arange(len(keep)), keep)), shape=(len(keep), D))
    proj = lambda A: (P @ A @ P.T).tocsr()
    dd = 2 if Ly > 1 else 1; M = m + dd*r
    H = sp.csr_matrix((len(keep), len(keep)), dtype=complex)
    for s in range(N):
        H += M*(proj(op({2*s: a.T@a})) - proj(op({2*s+1: a.T@a})))
    for l in range(nl):
        H += g*g/2*proj(op({nm+l: Esz@Esz}))
    for (b, rr, tt, ll) in plaq:
        Upl = proj(op({nm+b: Up}) @ op({nm+rr: Up}) @ op({nm+tt: Up.T}) @ op({nm+ll: Up.T}))
        H += -1/(2*g*g)*(Upl + Upl.T.conj())
    sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1])
    for l, (s, t_, d) in enumerate(links):
        A = r*sz - 1j*(sx if d == 0 else sy)
        U = op({nm+l: Up})
        for aa in range(2):
            for bb in range(2):
                if A[aa, bb] == 0: continue
                cd = op(c(2*s+aa)).T  # creation (real matrices)
                cb = op(c(2*t_+bb))
                T = proj(cd @ U @ cb) * (-0.5*A[aa, bb])
                H += T + T.T.conj()
    herm = abs(H - H.T.conj()).max()
    vals = np.sort(sla.eigsh(H, k=4, which='SA')[0]) if H.shape[0] > 50 else np.sort(np.linalg.eigvalsh(H.toarray()))[:4]
    # Loschmidt amplitude of bare Dirac sea with flux-free links (only if Q==0)
    vac_digits = [0, 1]*N + [1]*nl
    vi = 0
    for dg, d in zip(vac_digits, dims): vi = vi*d + dg
    loschmidt = None
    pos = np.searchsorted(keep, vi)
    if pos < len(keep) and keep[pos] == vi:
        v = np.zeros(len(keep), complex); v[pos] = 1
        w = sla.expm_multiply(-1j*t*H, v)
        loschmidt = w[pos]
    print(f"{Lx}x{Ly} {'per' if per else 'open'} Q={Qv} dim={H.shape[0]} herm={herm:.1e} E0..3={np.round(vals, 8).tolist()} L(t)={loschmidt}")

for args in [(2,1,False),(2,1,True),(3,1,True),(4,1,False),(2,2,False)]:
    run(*args)
run(3,1,False,Q=[1,0,-1])
run(2,2,True)
