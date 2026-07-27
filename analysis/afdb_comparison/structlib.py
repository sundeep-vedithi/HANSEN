"""
Self-contained mmCIF Ca parser + structural comparison metrics.

No external dependencies beyond numpy. Written for the HANSEN AFDB-vs-predictors
supplementary analysis.

Metrics implemented from primary definitions:
  TM-score      Zhang & Skolnick (2004) Proteins 57:702-710
  GDT-TS/GDT-HA Zemla (2003) Nucleic Acids Res 31:3370-3374
  lDDT          Mariani et al. (2013) Bioinformatics 29:2722-2728
  P-SEA SS      Labesse et al. (1997) CABIOS 13:291-295  (Ca-only assignment)
"""
import gzip
import math
import numpy as np

# ---------------------------------------------------------------- mmCIF parsing

THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V', 'MSE': 'M', 'SEC': 'U', 'PYL': 'O',
    'UNK': 'X',
}


def _tokenize_cif_line(line):
    """Whitespace split honouring single/double quotes used by mmCIF."""
    out, buf, quote = [], [], None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
                out.append(''.join(buf))
                buf = []
            else:
                buf.append(ch)
        elif ch in ("'", '"'):
            quote = ch
        elif ch.isspace():
            if buf:
                out.append(''.join(buf))
                buf = []
        else:
            buf.append(ch)
    if buf:
        out.append(''.join(buf))
    return out


def _open(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


def parse_ca(path, chain=None, model=1):
    """Parse the atom_site loop of an mmCIF file and return CA-level arrays.

    Returns dict with keys:
        resnum  (int64 [N])   auth_seq_id
        seq     (str  N)      one-letter sequence in file order
        xyz     (float64 [N,3])
        plddt   (float64 [N]) B_iso_or_equiv (pLDDT for all model sources here)
        chain   (str)         auth chain actually used
    Only the first model and a single chain are returned (all inputs here are
    monomers; if several chains exist the first encountered is used).
    """
    cols = {}
    order = []
    in_loop = False
    rows = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith('_atom_site.'):
                key = line.strip().split('.', 1)[1].split()[0]
                cols[key] = len(order)
                order.append(key)
                in_loop = True
                continue
            if in_loop:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    rows.append(line)
                elif line.startswith('#') or line.startswith('loop_') or line.startswith('_'):
                    if rows:
                        break
                    in_loop = bool(cols)
    if not rows:
        raise ValueError(f'no atom_site records parsed from {path}')

    def gi(tok, name, default=None):
        i = cols.get(name)
        if i is None or i >= len(tok):
            return default
        return tok[i]

    resnum, seq, xyz, bf = [], [], [], []
    used_chain = chain
    for line in rows:
        tok = _tokenize_cif_line(line)
        if gi(tok, 'group_PDB') not in ('ATOM', 'HETATM'):
            continue
        if gi(tok, 'label_atom_id') != 'CA':
            continue
        if gi(tok, 'type_symbol') != 'C':
            continue
        mdl = gi(tok, 'pdbx_PDB_model_num', '1')
        if mdl not in ('.', '?', None) and int(mdl) != model:
            continue
        alt = gi(tok, 'label_alt_id', '.')
        if alt not in ('.', '?', 'A', None):
            continue
        ch = gi(tok, 'auth_asym_id') or gi(tok, 'label_asym_id')
        if used_chain is None:
            used_chain = ch
        if ch != used_chain:
            continue
        comp = gi(tok, 'auth_comp_id') or gi(tok, 'label_comp_id')
        one = THREE2ONE.get(comp, 'X')
        rn = gi(tok, 'auth_seq_id') or gi(tok, 'label_seq_id')
        try:
            rn = int(rn)
        except (TypeError, ValueError):
            continue
        try:
            x = float(gi(tok, 'Cartn_x'))
            y = float(gi(tok, 'Cartn_y'))
            z = float(gi(tok, 'Cartn_z'))
        except (TypeError, ValueError):
            continue
        try:
            b = float(gi(tok, 'B_iso_or_equiv'))
        except (TypeError, ValueError):
            b = float('nan')
        resnum.append(rn)
        seq.append(one)
        xyz.append((x, y, z))
        bf.append(b)

    if not resnum:
        raise ValueError(f'no CA atoms parsed from {path}')
    return {
        'resnum': np.asarray(resnum, dtype=np.int64),
        'seq': ''.join(seq),
        'xyz': np.asarray(xyz, dtype=np.float64),
        'plddt': np.asarray(bf, dtype=np.float64),
        'chain': used_chain,
    }


# ------------------------------------------------------------------- alignment

def align_by_resnum(a, b):
    """Residue correspondence by auth_seq_id (identical UniProt numbering)."""
    ib = {int(r): i for i, r in enumerate(b['resnum'])}
    ia_idx, ib_idx = [], []
    for i, r in enumerate(a['resnum']):
        j = ib.get(int(r))
        if j is not None:
            ia_idx.append(i)
            ib_idx.append(j)
    return np.asarray(ia_idx, dtype=np.int64), np.asarray(ib_idx, dtype=np.int64)


def needleman_wunsch(s1, s2, match=1.0, mismatch=-1.0, gap=-2.0):
    """Global alignment; returns index pairs of aligned (non-gap) positions."""
    n, m = len(s1), len(s2)
    score = np.zeros((n + 1, m + 1), dtype=np.float32)
    score[:, 0] = np.arange(n + 1) * gap
    score[0, :] = np.arange(m + 1) * gap
    ptr = np.zeros((n + 1, m + 1), dtype=np.int8)   # 0 diag, 1 up, 2 left
    ptr[:, 0] = 1
    ptr[0, :] = 2
    a1 = np.frombuffer(s1.encode(), dtype=np.uint8)
    a2 = np.frombuffer(s2.encode(), dtype=np.uint8)
    for i in range(1, n + 1):
        sim = np.where(a1[i - 1] == a2, match, mismatch)
        prev = score[i - 1]
        row = score[i]
        prow = ptr[i]
        best_left = -np.inf
        for j in range(1, m + 1):
            d = prev[j - 1] + sim[j - 1]
            u = prev[j] + gap
            left = row[j - 1] + gap
            if d >= u and d >= left:
                row[j] = d
                prow[j] = 0
            elif u >= left:
                row[j] = u
                prow[j] = 1
            else:
                row[j] = left
                prow[j] = 2
        del best_left
    i, j = n, m
    p1, p2 = [], []
    while i > 0 or j > 0:
        d = ptr[i, j]
        if d == 0:
            i -= 1
            j -= 1
            p1.append(i)
            p2.append(j)
        elif d == 1:
            i -= 1
        else:
            j -= 1
    return np.asarray(p1[::-1], dtype=np.int64), np.asarray(p2[::-1], dtype=np.int64)


# ------------------------------------------------------------- superposition

def kabsch_rt(P, Q):
    """Rotation+translation minimising RMSD of P onto Q. Returns (R, t)."""
    pc = P.mean(axis=0)
    qc = Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = qc - R @ pc
    return R, t


def apply_rt(P, R, t):
    return (R @ P.T).T + t


def rmsd_after_fit(P, Q):
    R, t = kabsch_rt(P, Q)
    d = apply_rt(P, R, t) - Q
    return float(np.sqrt((d ** 2).sum() / len(P)))


# ------------------------------------------------------------------ TM-score

def _d0(L):
    if L <= 15:
        return 0.5
    return max(0.5, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)


def tm_score(P, Q, L_norm=None, n_iter=20):
    """Sequence-dependent TM-score of model P against reference Q.

    P, Q are [N,3] arrays of already-corresponded CA coordinates.
    Normalised by L_norm (default: N).  Implements the Zhang-Skolnick
    iterative fragment-seed search.
    Returns (tm, rmsd_of_best_superposed_subset, R, t).
    """
    N = len(P)
    if N < 3:
        return 0.0, float('nan'), np.eye(3), np.zeros(3)
    L = int(L_norm) if L_norm else N
    d0 = _d0(L)
    d0sq = d0 * d0

    best_tm = -1.0
    best_R, best_t = np.eye(3), np.zeros(3)

    # seed fragment lengths, classic schedule
    seeds = []
    Lf = N
    while Lf >= 4:
        seeds.append(Lf)
        Lf //= 2
    if seeds[-1] != 4:
        seeds.append(4)

    for Lf in seeds:
        stride = max(1, Lf // 2) if Lf < N else 1
        for start in range(0, N - Lf + 1, stride):
            idx = np.arange(start, start + Lf)
            for _ in range(n_iter):
                if len(idx) < 3:
                    break
                R, t = kabsch_rt(P[idx], Q[idx])
                Pt = apply_rt(P, R, t)
                dsq = ((Pt - Q) ** 2).sum(axis=1)
                tm = float((1.0 / (1.0 + dsq / d0sq)).sum() / L)
                if tm > best_tm:
                    best_tm, best_R, best_t = tm, R, t
                # re-select residues within a growing cutoff
                dcut = max(d0, 3.0)
                new = np.where(dsq < dcut * dcut)[0]
                grow = dcut
                while len(new) < 3 and grow < 20.0:
                    grow += 0.5
                    new = np.where(dsq < grow * grow)[0]
                if len(new) == len(idx) and np.array_equal(new, idx):
                    break
                idx = new
    return best_tm, best_R, best_t


# --------------------------------------------------------------------- GDT

def gdt(P, Q, cutoffs=(1.0, 2.0, 4.0, 8.0), n_iter=30):
    """GDT score: mean over cutoffs of the largest fraction of residues that can
    be superposed within that cutoff (iterative superposition per cutoff)."""
    N = len(P)
    fracs = []
    for c in cutoffs:
        best = 0
        for Lf in (N, max(3, N // 2), max(3, N // 4), 7, 3):
            for start in (0, max(0, (N - Lf) // 2), max(0, N - Lf)):
                idx = np.arange(start, min(N, start + Lf))
                if len(idx) < 3:
                    continue
                for _ in range(n_iter):
                    R, t = kabsch_rt(P[idx], Q[idx])
                    d = np.linalg.norm(apply_rt(P, R, t) - Q, axis=1)
                    new = np.where(d < c)[0]
                    if len(new) > best:
                        best = len(new)
                    if len(new) < 3 or np.array_equal(new, idx):
                        break
                    idx = new
        fracs.append(best / N)
    return float(np.mean(fracs)), fracs


# -------------------------------------------------------------------- lDDT

def lddt_ca(P, Q, inclusion_radius=15.0, thresholds=(0.5, 1.0, 2.0, 4.0),
            sep=0):
    """Superposition-free local Distance Difference Test on CA atoms.

    P = model coords, Q = reference coords, already corresponded.
    Returns (global_lddt, per_residue_lddt [N]).
    """
    N = len(P)
    dQ = np.linalg.norm(Q[:, None, :] - Q[None, :, :], axis=-1)
    dP = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    mask = (dQ < inclusion_radius)
    np.fill_diagonal(mask, False)
    if sep:
        ii = np.arange(N)
        mask &= (np.abs(ii[:, None] - ii[None, :]) > sep)
    diff = np.abs(dP - dQ)
    preserved = np.zeros_like(diff, dtype=np.float64)
    for th in thresholds:
        preserved += (diff < th)
    preserved /= len(thresholds)
    denom = mask.sum(axis=1)
    per_res = np.where(denom > 0, (preserved * mask).sum(axis=1) / np.maximum(denom, 1), np.nan)
    tot = mask.sum()
    glob = float((preserved * mask).sum() / tot) if tot else float('nan')
    return glob, per_res


# ------------------------------------------------------- secondary structure

def psea_ss(xyz):
    """Ca-only secondary structure assignment (P-SEA; Labesse et al. 1997).

    Returns a string of 'H' (helix), 'E' (strand), 'C' (coil), one per CA.

    Reference values from the paper (mean +/- tolerance):
      helix  d(i,i+2) 5.5+/-0.5   d(i,i+3) 5.3+/-0.5   d(i,i+4) 6.4+/-0.6
             alpha 89.4+/-12 deg  tau 49.9+/-7.6 deg
      strand d(i,i+2) 6.7+/-0.6   d(i,i+3) 9.9+/-0.9   d(i,i+4) 12.4+/-1.1
             alpha 124.7+/-14.5 deg  tau -170.0+/-45 deg
    An element is assigned when either the distance criteria or the angle
    criteria are satisfied, as in the original three-distance/two-angle scheme.
    """
    N = len(xyz)
    ss = np.array(['C'] * N)
    if N < 6:
        return ''.join(ss)

    def d(i, j):
        if j >= N or i < 0:
            return np.inf
        return float(np.linalg.norm(xyz[i] - xyz[j]))

    def bond_angle(i):
        if i + 2 >= N:
            return np.nan
        v1 = xyz[i] - xyz[i + 1]
        v2 = xyz[i + 2] - xyz[i + 1]
        c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
        return math.degrees(math.acos(max(-1.0, min(1.0, c))))

    def torsion(i):
        if i + 3 >= N:
            return np.nan
        b1 = xyz[i + 1] - xyz[i]
        b2 = xyz[i + 2] - xyz[i + 1]
        b3 = xyz[i + 3] - xyz[i + 2]
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        m = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-9))
        return math.degrees(math.atan2(float(np.dot(m, n2)), float(np.dot(n1, n2))))

    # --- helices
    for i in range(N - 4):
        d2, d3, d4 = d(i, i + 2), d(i, i + 3), d(i, i + 4)
        dist_ok = (5.0 <= d2 <= 6.0) and (4.8 <= d3 <= 5.8) and (5.8 <= d4 <= 7.0)
        a, t = bond_angle(i), torsion(i)
        ang_ok = (a == a and t == t and 77.4 <= a <= 101.4 and 42.3 <= t <= 57.5)
        if dist_ok or ang_ok:
            ss[i:i + 4] = 'H'
    # --- strands
    for i in range(N - 4):
        d2, d3, d4 = d(i, i + 2), d(i, i + 3), d(i, i + 4)
        dist_ok = (6.1 <= d2 <= 7.3) and (9.0 <= d3 <= 10.8) and (11.3 <= d4 <= 13.5)
        a, t = bond_angle(i), torsion(i)
        ang_ok = (a == a and t == t and 110.2 <= a <= 139.2 and
                  (t <= -125.0 or t >= 145.0))
        if dist_ok or ang_ok:
            for k in range(i, min(N, i + 3)):
                if ss[k] != 'H':
                    ss[k] = 'E'
    # --- drop isolated single-residue elements (P-SEA post-filter)
    out = ss.copy()
    for i in range(N):
        if ss[i] != 'C':
            left = ss[i - 1] if i > 0 else 'C'
            right = ss[i + 1] if i < N - 1 else 'C'
            if left != ss[i] and right != ss[i]:
                out[i] = 'C'
    return ''.join(out)


def radius_of_gyration(xyz):
    c = xyz.mean(axis=0)
    return float(np.sqrt(((xyz - c) ** 2).sum(axis=1).mean()))


def pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float('nan')
    a, b = a[m], b[m]
    if a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float('nan')
    ra = np.argsort(np.argsort(a[m])).astype(float)
    rb = np.argsort(np.argsort(b[m])).astype(float)
    return pearson(ra, rb)
