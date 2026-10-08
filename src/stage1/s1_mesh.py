"""
Structured hexahedral mesh for the Stage 1 single-pad unit cell (Ayoub 2022 reproduction).

Geometry (uMKS, um), z = 0 at the Si / SiO2 interface:
    Si slab        -T_SI .. 0           (substrate equivalent, see spec B4)
    SiO2 under pad  0    .. Z_PB        (0.30 um)
    pad layer       Z_PB .. Z_PT        (0.85 um; Cu column in SiO2)
    SiN cap         Z_PT .. Z_PT+0.06   (optional)
Lateral cell 2.9 x 2.9 um, pad 0.3 x 0.3 um centred (square) or D = 0.3 um (circular O-grid).

The mesh is generated here (not by MAPDL meshing) so that node IDs of the
probe points, the sidewall node set and the duplicated sidewall nodes of the
sliding variant are known exactly and identically in PyMAPDL and batch runs.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# ------------------------------------------------------------------ geometry
L_CELL = 2.9          # um  (2.6 um edge-to-edge spacing + 0.3 um pad)
A_PAD = 0.15          # um  half width of the pad (square) / radius (circle)
C_PAD = L_CELL / 2.0
Z_PB = 0.30           # um  pad bottom (SiO2 below pad: ASSUMPTION A3)
H_PAD = 0.85          # um  pad height (Fig. 1a)
Z_PT = Z_PB + H_PAD   # 1.15 um
T_CAP = 0.06          # um  SiN cap (paper Sec. 2)
T_SI = 2.0            # um  Si slab thickness (substrate equivalent)

MAT_CU, MAT_OX, MAT_SI, MAT_SIN = 1, 2, 3, 4

# mesh-refinement levels (G2): in-pad element size and grading limits
LEVELS = {
    "L0": dict(n_pad=6,  growth=1.30, h_max_lat=0.30, h_max_z=0.30),
    "L1": dict(n_pad=8,  growth=1.28, h_max_lat=0.26, h_max_z=0.26),
    "L2": dict(n_pad=10, growth=1.25, h_max_lat=0.22, h_max_z=0.22),
    "L3": dict(n_pad=12, growth=1.22, h_max_lat=0.20, h_max_z=0.20),
}
STUDENT_LIMIT = 128_000


def _graded(length, h0, growth, h_max):
    """Increments starting at h0, growing geometrically, capped, scaled to sum to length."""
    inc = []
    h = h0
    while sum(inc) < length - 1e-12:
        inc.append(min(h, h_max))
        h *= growth
    inc = np.array(inc)
    # merge a tiny last increment into its neighbour, then scale to exact length
    if len(inc) > 1 and inc[-1] < 0.5 * inc[-2]:
        inc = inc[:-1]
    return inc * (length / inc.sum())


def _symmetric_graded(length, h_end, growth, h_max):
    """Fine at both ends, coarse in the middle."""
    half = _graded(length / 2.0, h_end, growth, h_max)
    return np.concatenate([half, half[::-1]])


def lateral_axis(n_pad, growth, h_max):
    h_p = 2 * A_PAD / n_pad
    outside = _graded(C_PAD - A_PAD, h_p, growth, h_max)        # from pad edge outwards
    inc = np.concatenate([outside[::-1], np.full(n_pad, h_p), outside])
    x = np.concatenate([[0.0], np.cumsum(inc)])
    x[-1] = L_CELL
    i0 = len(outside)                                           # first pad grid line
    return x, i0, i0 + n_pad


def vertical_axis(n_pad, growth, h_max, cap=True, t_si=T_SI):
    h_p = 2 * A_PAD / n_pad
    si = _graded(t_si, h_p, growth, h_max)[::-1]                 # coarse at bottom
    ox = _symmetric_graded(Z_PB, h_p, growth, h_max)
    pad = _symmetric_graded(H_PAD, h_p, growth, 2.5 * h_p)       # fine at top/bottom
    segs = [si, ox, pad]
    if cap:
        n_cap = max(2, int(math.ceil(T_CAP / h_p)))
        segs.append(np.full(n_cap, T_CAP / n_cap))
    inc = np.concatenate(segs)
    z = np.concatenate([[-t_si], -t_si + np.cumsum(inc)])
    k_si = len(si)
    k_pb = k_si + len(ox)
    k_pt = k_pb + len(pad)
    # snap interfaces exactly
    z[k_si], z[k_pb], z[k_pt] = 0.0, Z_PB, Z_PT
    if cap:
        z[-1] = Z_PT + T_CAP
    return z, k_si, k_pb, k_pt


def _ogrid_map(x, y, core=0.5, s_out=3.0):
    """Square-to-circle O-grid blend about the pad centre (normalised by A_PAD).

    s = max(|u|,|v|) (square radius). s <= core: identity; core..1: square->circle;
    1..s_out: circle->square; s >= s_out: identity. Pad boundary s = 1 maps to the
    circle of radius A_PAD (diameter 300 nm)."""
    u, v = (x - C_PAD) / A_PAD, (y - C_PAD) / A_PAD
    s = np.maximum(np.abs(u), np.abs(v))
    with np.errstate(invalid="ignore", divide="ignore"):
        su, sv = np.where(s > 0, u / s, 0.0), np.where(s > 0, v / s, 0.0)   # unit-square dir
        r = np.hypot(su, sv)
        cu, cv = np.where(r > 0, su / r, 0.0), np.where(r > 0, sv / r, 0.0) # unit-circle dir
    t_in = np.clip((s - core) / (1.0 - core), 0.0, 1.0)
    t_out = np.clip((s - 1.0) / (s_out - 1.0), 0.0, 1.0)
    # weight of the circular direction: 0 at core, 1 at s = 1, back to 0 at s_out
    w = np.where(s <= 1.0, t_in, 1.0 - t_out)
    du = s * ((1 - w) * su + w * cu)
    dv = s * ((1 - w) * sv + w * cv)
    return C_PAD + A_PAD * du, C_PAD + A_PAD * dv


@dataclass
class Mesh:
    level: str
    shape: str
    cap: bool
    sliding: bool
    x: np.ndarray = field(repr=False)
    z: np.ndarray = field(repr=False)
    nodes: np.ndarray = field(repr=False)          # (N, 4): id, X, Y, Z
    elems: dict = field(repr=False)                # mat -> (M, 8) node ids
    ce_pairs: list = field(repr=False)             # (cu_node, ox_node, nx, ny)
    probes: dict = field(default_factory=dict)
    sidewall: list = field(default_factory=list, repr=False)   # (node, nx, ny) Cu side
    cu_top: list = field(default_factory=list, repr=False)
    n_nodes: int = 0
    n_elems: int = 0


def build_mesh(level="L2", shape="square", cap=True, sliding=False, t_si=T_SI):
    p = LEVELS[level]
    x, i0, i1 = lateral_axis(p["n_pad"], p["growth"], p["h_max_lat"])
    z, k_si, k_pb, k_pt = vertical_axis(p["n_pad"], p["growth"], p["h_max_z"], cap, t_si)
    nx = ny = len(x) - 1
    nz = len(z) - 1
    X, Y = np.meshgrid(x, x, indexing="ij")
    if shape == "circle":
        X, Y = _ogrid_map(X, Y)
    elif shape != "square":
        raise ValueError(shape)

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    I, J, K = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1), np.arange(nz + 1), indexing="ij")
    ids = nid(I, J, K)
    nodes = np.column_stack([ids.ravel(), X[I, J].ravel(), Y[I, J].ravel(), z[K].ravel()])
    nodes = nodes[np.argsort(nodes[:, 0])]

    # elements
    ei, ej, ek = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    ei, ej, ek = ei.ravel(), ej.ravel(), ek.ravel()
    conn = np.column_stack([nid(ei, ej, ek), nid(ei + 1, ej, ek), nid(ei + 1, ej + 1, ek), nid(ei, ej + 1, ek),
                            nid(ei, ej, ek + 1), nid(ei + 1, ej, ek + 1), nid(ei + 1, ej + 1, ek + 1), nid(ei, ej + 1, ek + 1)])
    in_pad_xy = (ei >= i0) & (ei < i1) & (ej >= i0) & (ej < i1)
    mat = np.where(ek < k_si, MAT_SI,
          np.where(ek < k_pb, MAT_OX,
          np.where(ek < k_pt, np.where(in_pad_xy, MAT_CU, MAT_OX), MAT_SIN)))

    # sidewall nodes: pad boundary lines, strictly between pad bottom and top
    sw = []
    for k in range(k_pb, k_pt + 1):
        for i in range(i0, i1 + 1):
            for j in range(i0, i1 + 1):
                on_x = i in (i0, i1)
                on_y = j in (i0, i1)
                if not (on_x or on_y):
                    continue
                if shape == "square":
                    nxv = (-1.0 if i == i0 else 1.0) if on_x else 0.0
                    nyv = (-1.0 if j == i0 else 1.0) if on_y else 0.0
                else:
                    dx, dy = X[i, j] - C_PAD, Y[i, j] - C_PAD
                    rr = math.hypot(dx, dy)
                    nxv, nyv = dx / rr, dy / rr
                nn = math.hypot(nxv, nyv)
                sw.append((nid(i, j, k), k, i, j, nxv / nn, nyv / nn, on_x and on_y))

    ce_pairs = []
    sidewall = [(n, nxv, nyv) for (n, k, i, j, nxv, nyv, corner) in sw]
    if sliding:
        # Frictionless, no-separation sidewall. Cu-side copies of sidewall nodes are created and
        # tied to the SiO2 node by the NORMAL displacement only (UZ and tangential free):
        #   - sidewall face nodes (strictly between pad bottom and top): tie the face normal;
        #   - vertical corner-edge nodes (square pad): two faces meet -> tie UX and UY, UZ free;
        #   - top rim (k = k_pt): duplicated ONLY when there is no cap (with a bonded cap the cap
        #     joins Cu and SiO2 at the rim, so the shared node is physical); corners tie UX and UY;
        #   - bottom rim (k = k_pb) stays shared: the pad bottom is bonded by design.
        dup = {}
        next_id = int(nodes[:, 0].max()) + 1
        new_rows = []
        for (n, k, i, j, nxv, nyv, corner) in sw:
            if k == k_pb:
                continue
            if k == k_pt and cap:
                continue
            dup[n] = next_id
            row = nodes[n - 1].copy(); row[0] = next_id
            new_rows.append(row)
            if shape == "square" and corner:
                ce_pairs.append((next_id, n, 1.0, 0.0))
                ce_pairs.append((next_id, n, 0.0, 1.0))
            else:
                ce_pairs.append((next_id, n, nxv, nyv))
            next_id += 1
        nodes = np.vstack([nodes, np.array(new_rows)])
        cu = mat == MAT_CU
        sub = conn[cu]
        remap = np.vectorize(lambda v: dup.get(int(v), int(v)))
        conn[cu] = remap(sub)
        sidewall = [(dup.get(n, n), nxv, nyv) for (n, nxv, nyv) in sidewall]
    else:
        dup = {}

    elems = {m: conn[mat == m] for m in (MAT_CU, MAT_OX, MAT_SI, MAT_SIN) if np.any(mat == m)}

    ic = (i0 + i1) // 2
    probes = dict(
        cu_top_center=nid(ic, ic, k_pt),
        cu_top_edge=nid(i1, ic, k_pt),          # mid-edge on the +x sidewall (square)
        cu_bottom_center=nid(ic, ic, k_pb),
        si_top_corner=nid(0, 0, k_si),
        surface_corner=nid(0, 0, nz),           # cap top (or SiO2 top if no cap), far field
        surface_center=nid(ic, ic, nz),
        surface_edge_mid=nid(0, ic, nz),     # top surface (cap/dielectric) at the cell edge, mid-side
    )
    cu_top = [dup.get(nid(i, j, k_pt), nid(i, j, k_pt)) for i in range(i0, i1 + 1) for j in range(i0, i1 + 1)]
    for key in ("cu_top_center", "cu_top_edge", "cu_bottom_center"):
        probes[key] = dup.get(probes[key], probes[key])
    m = Mesh(level, shape, cap, sliding, x, z, nodes, elems, ce_pairs, probes, sidewall, cu_top)
    m.n_nodes = len(nodes)
    m.n_elems = int(sum(len(v) for v in elems.values()))
    m.k = dict(k_si=k_si, k_pb=k_pb, k_pt=k_pt, nz=nz, i0=i0, i1=i1, nx=nx)
    m.t_si = t_si
    return m


def element_quality(m):
    """Min/max corner-angle proxy for the lateral (x-y) faces of all elements."""
    xy = {int(r[0]): (r[1], r[2]) for r in m.nodes}
    worst = 0.0
    for conn in m.elems.values():
        for e in conn[:: max(1, len(conn) // 4000)]:
            q = [xy[int(v)] for v in e[:4]]
            for a in range(4):
                p0, p1, p2 = np.array(q[a - 1]), np.array(q[a]), np.array(q[(a + 1) % 4])
                v1, v2 = p0 - p1, p2 - p1
                ang = math.degrees(math.acos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)))
                worst = max(worst, ang)
    return worst
