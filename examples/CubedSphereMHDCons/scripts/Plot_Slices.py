#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot combined Tecplot ASCII slice (BLOCK packing) with multiple ZONEs
into a single domain using cell centers (for cell-centered variables).

Assumptions:
- VARIABLES = ["X","Y","Z", physics...]
- First 3 (X,Y,Z) are node-centered (size I*J each)
- Variables 4..N are cell-centered (size (I-1)*(J-1) each)
- ZONE ... DATAPACKING=BLOCK, VARLOCATION=([4-12]=CELLCENTERED)
- Files may have many ZONEs (e.g., AMR blocks)

Compatible with AMR: concatenates zones and prunes overlaps by keeping the finest level.
"""

import os
import re
from typing import List, Dict, Any, Tuple
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.patches import Circle
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# --- physics conversion params (edit as needed) ---
AU_CM       = 1.495978707e13  # AU in cm
START_ISO   = "2024-03-16T12:00:00"  # simulation t=0 (UTC)

#--- address variables ---
# slice_folder = "/Users/talwindersingh/Desktop/My_Computer/Work/UAH/Current_projects/Ron/With_HelioCubed_single_frame"
slice_folder = "/Users/talwindersingh/Desktop/My_Computer/Work/UAH/Current_projects/Ron/With_HelioCubed"

out_dir = slice_folder + "/Images_inertial"

#-- Corotation needed? ---
COROTATE = False

#-- process all or new only ---
PROCESS_ALL = False

# -------------------------------
# Radial scaling powers per variable
# -------------------------------
RADIAL_SCALE = {
    "Vr": 0,
    "Vt": 0,
    "Vp": 0,
    "Br": 2,
    "Bt": 2,
    "Bp": 2,
    "density": 2,
    "T": 1,
    "P": 2.1
}

UNIT_SCALE = {
    "Br": 0.1e6,  # G -> nT
    "Bt": 0.1e6,  # G -> nT
    "Bp": 0.1e6,  # G -> nT
    "Vr": 1e-5, # cm/s -> km/s
    "Vt": 1e-5, # cm/s -> km/s
    "Vp": 1e-5, # cm/s -> km/s
    "density": 1.0, # g/cm³
    "T": 1.0,
    "P": 1.0e12,  # dyn/cm² -> pico dyn/cm²
}

LABELS = {
    "Vr": "Vr [km/s]", "Vp": "Vp [km/s]", "Vt": "Vt [km/s]",
    "Br": "Br_scaled [nT]",   "Bp": "Bp_scaled [nT]",   "Bt": "Bt_scaled [nT]",
    "density": "ρ_scaled [g/cm³]", "T": "T [K]", "P": "p_scaled [pico dyn/cm²]"
}

# -------------------------------
# Colormaps and vmin/vmax settings
# -------------------------------
COLOR_SETTINGS = {
    # Velocities
    "Vr": {"cmap": "rainbow", "vmin": 200.0, "vmax": 700.0},
    "Vt": {"cmap": "rainbow", "vmin": -10.0, "vmax": 10.0},
    "Vp": {"cmap": "rainbow", "vmin": -10.0, "vmax": 10.0},

    # Magnetic Fields
    "Br": {"cmap": "bwr", "vmin": -3.0, "vmax": 3.0},
    "Bp": {"cmap": "bwr", "vmin": -1.0, "vmax": 1.0},
    "Bt": {"cmap": "bwr", "vmin": -0.2, "vmax": 0.2},

    # Mass density
    "density": {"cmap": "terrain_r", "vmin": 0.0, "vmax": 0.25e-22},

    # Temperature
    "T": {"cmap": "hot", "vmin": 0.0, "vmax": 1.5e5},

    # Pressure
    "P": {"cmap": "plasma", "vmin": 0.0, "vmax": 800},
}

def cm_to_AU(x, pos):
    return f"{x / AU_CM:.2f}"

# =========================
# Geometry cache + builders
# =========================
GEOM_CACHE = {}  # key = (filename, reshape_order, corotate, rot_rate)

def _reshape_conn_to_e4(conn: np.ndarray) -> np.ndarray:
    """Ensure FE connectivity shaped (E,4); accept (E,4), (4,E), or flat length 4E."""
    conn = np.asarray(conn)
    if conn.ndim == 1:
        if conn.size % 4 != 0:
            raise ValueError(f"Flat connectivity length {conn.size} not divisible by 4")
        return conn.reshape((-1, 4))
    if conn.ndim == 2:
        if conn.shape[1] == 4:
            return conn
        if conn.shape[0] == 4:
            return conn.T
    raise ValueError(f"Unsupported connectivity shape {conn.shape}")

def _normalize_conn_indices(conn: np.ndarray, n_nodes: int) -> np.ndarray:
    """Return 0-based (E,4) connectivity; handles 1-based Tecplot indices."""
    conn = _reshape_conn_to_e4(conn).astype(int, copy=False)
    cmin, cmax = conn.min(), conn.max()
    if cmin >= 1 and cmax <= n_nodes:
        conn0 = conn - 1
    elif cmin >= 0 and cmax <= n_nodes - 1:
        conn0 = conn
    else:
        raise ValueError(f"Connectivity out of range: min={cmin}, max={cmax}, n_nodes={n_nodes}")
    if (conn0 < 0).any() or (conn0 >= n_nodes).any():
        bad = conn0[(conn0 < 0) | (conn0 >= n_nodes)][:8]
        raise IndexError(f"Connectivity out of bounds after normalization, e.g. {bad}")
    return conn0

def _fe_nodes(z):
    """Get (Xn, Yn) from FE zone that stores nodes as dict or tuple/list."""
    nd = z.get('nodes')
    if isinstance(nd, dict):
        Xn = np.asarray(nd.get('X', nd.get('x')))
        Yn = np.asarray(nd.get('Y', nd.get('y')))
        if Xn is None or Yn is None:
            raise KeyError("FE nodes dict missing 'X'/'Y'")
        return Xn, Yn
    if isinstance(nd, (list, tuple)) and len(nd) >= 2:
        return np.asarray(nd[0]), np.asarray(nd[1])
    raise TypeError(f"Unsupported FE nodes type: {type(nd)}")

def _rotate_xy(x: np.ndarray, y: np.ndarray, theta: float) -> Tuple[np.ndarray, np.ndarray]:
    if theta == 0.0:
        return x, y
    c, s = np.cos(theta), np.sin(theta)
    return c*x - s*y, s*x + c*y

def _rotation_angle_since_start(phys_dt: datetime, rot_rate_deg_per_day: float, start_iso: str) -> float:
    """Angle sun has rotated since start, in radians (CCW positive)."""
    t0 = datetime.fromisoformat(start_iso)
    dt_days = (phys_dt - t0).total_seconds() / 86400.0
    return np.deg2rad(rot_rate_deg_per_day * dt_days)

def prepare_geometry(
    filename: str,
    reshape_order: str = "C",
    corotate: bool = True,
    rot_rate_deg_per_day: float = 14.1843971631,
    inner_radius: Optional[float] = None,
    outer_radius: Optional[float] = None,
):
    """
    Read the slice once; compute:
      - cell centers x/y (combined across zones)
      - cell radii Rc (from node radii; smoother)
      - triangulation + annulus mask
      - time/step + rotation angle used
    Return a dict 'geom' reused for all variables.
    """
    ds = read_tecplot_slice_ascii(filename)
    zones = ds['zones']
    variables = ds['variables']

    # time + rotation angle (corotating frame)
    step = -1; t_code = None; phys_dt = None; theta_corotate = 0.0
    try:
        t_code, step = parse_title_time_step(filename)
        phys_dt = code_time_to_physical(t_code)
        theta_inertial = _rotation_angle_since_start(phys_dt, rot_rate_deg_per_day, START_ISO)
        theta_corotate = -theta_inertial if corotate else 0.0
    except Exception:
        pass

    node_x_all, node_y_all = [], []
    xc_all, yc_all, rc_all = [], [], []

    for z in zones:
        props = z.get('props', {})
        zonetype = (props.get('ZONETYPE', '') or '').upper()

        # Structured
        if ('I' in z and 'J' in z) or ('I' in props and 'J' in props):
            I = int(z.get('I', props.get('I')))
            J = int(z.get('J', props.get('J')))

            X = np.reshape(z['blocks']['node_blocks'][0], (J, I), order=reshape_order)
            Y = np.reshape(z['blocks']['node_blocks'][1], (J, I), order=reshape_order)
            X, Y = _rotate_xy(X, Y, theta_corotate)

            node_x_all.append(X.ravel()); node_y_all.append(Y.ravel())

            XC = 0.25*(X[:-1,:-1] + X[1:,:-1] + X[:-1,1:] + X[1:,1:])
            YC = 0.25*(Y[:-1,:-1] + Y[1:,:-1] + Y[:-1,1:] + Y[1:,1:])
            Rn = np.hypot(X, Y)
            RC = 0.25*(Rn[:-1,:-1] + Rn[1:,:-1] + Rn[:-1,1:] + Rn[1:,1:])

            xc_all.append(XC.ravel()); yc_all.append(YC.ravel()); rc_all.append(RC.ravel())

        # FEQUADRILATERAL
        elif zonetype == 'FEQUADRILATERAL' or ('conn' in z):
            Xn, Yn = _fe_nodes(z)
            Xn, Yn = _rotate_xy(Xn, Yn, theta_corotate)
            node_x_all.append(Xn.ravel()); node_y_all.append(Yn.ravel())

            conn0 = _normalize_conn_indices(np.asarray(z['conn']), n_nodes=Xn.size)
            Xc = Xn[conn0].mean(axis=1)
            Yc = Yn[conn0].mean(axis=1)
            Rn = np.hypot(Xn, Yn)
            RC = Rn[conn0].mean(axis=1)

            xc_all.append(Xc.ravel()); yc_all.append(Yc.ravel()); rc_all.append(RC.ravel())

        else:
            raise KeyError("Unknown zone type (need I,J or FEQUADRILATERAL).")

    node_x = np.concatenate(node_x_all); node_y = np.concatenate(node_y_all)
    x_all  = np.concatenate(xc_all);     y_all  = np.concatenate(yc_all)
    rc_all = np.concatenate(rc_all)

    node_r = np.hypot(node_x, node_y)
    if inner_radius is None: inner_radius = float(node_r.min())
    if outer_radius is None: outer_radius = float(node_r.max())

    # triangulate + annulus mask once
    triang = tri.Triangulation(x_all, y_all)
    r_vert = np.hypot(x_all, y_all)
    tri_ix = triang.triangles
    mask = (r_vert[tri_ix] < inner_radius).any(axis=1) | (r_vert[tri_ix] > outer_radius).any(axis=1)
    triang.set_mask(mask)

    return {
        "filename": filename,
        "variables": variables,
        "x": x_all, "y": y_all, "rc": rc_all,
        "inner": inner_radius, "outer": outer_radius,
        "tri": triang,
        "theta": theta_corotate,
        "step": step, "t_code": t_code, "phys_dt": phys_dt,
        # keep zones to fetch cell data later
        "zones": zones,
        "reshape_order": reshape_order,
    }

def get_values_for_var(geom: dict, varname: str, cell_order: str = "C", rotate_components: bool = False) -> np.ndarray:
    """
    Pull cell-centered values for 'varname' across zones in the SAME order
    as geom['x']/['y']/['rc'] were assembled.
    """
    vals = []

    def rot_xy(a, b, th):
        if th == 0.0: return a, b
        c, s = np.cos(th), np.sin(th); return c*a - s*b, s*a + c*b

    want_vec = rotate_components and varname in ("Vx", "Vy", "Bx", "By")
    family = "V" if varname.startswith("V") else ("B" if varname.startswith("B") else None)
    pair = (family + ("y" if varname.endswith("x") else "x")) if (want_vec and family) else None

    reshape_order = geom["reshape_order"]

    for z in geom["zones"]:
        props = z.get('props', {})
        zonetype = (props.get('ZONETYPE', '') or '').upper()

        if ('I' in z and 'J' in z) or ('I' in props and 'J' in props):  # structured
            I = int(z.get('I', props.get('I')))
            J = int(z.get('J', props.get('J')))
            V = np.reshape(z['blocks']['cell_blocks'][varname], (J-1, I-1), order=cell_order)
            if want_vec and (pair in z['blocks']['cell_blocks']):
                W = np.reshape(z['blocks']['cell_blocks'][pair], (J-1, I-1), order=cell_order)
                if varname.endswith("x"):
                    V, _ = rot_xy(V, W, geom["theta"])
                else:
                    _, V = rot_xy(W, V, geom["theta"])
            vals.append(V.ravel())
        elif zonetype == 'FEQUADRILATERAL' or ('conn' in z):
            cells = z.get('cells', {})
            if varname not in cells:
                raise KeyError(f"FE cells missing '{varname}'")
            V = np.asarray(cells[varname]).ravel()
            if want_vec and (pair in cells):
                W = np.asarray(cells[pair]).ravel()
                if varname.endswith("x"):
                    V, _ = rot_xy(V, W, geom["theta"])
                else:
                    _, V = rot_xy(W, V, geom["theta"])
            vals.append(V)
        else:
            raise KeyError("Unknown zone type while fetching values.")

    return np.concatenate(vals)








# TITLE = "t = 1.23456789E-03 step = 12345"
TITLE_RE = re.compile(
    r'TITLE\s*=\s*"\s*t\s*=\s*([0-9Ee+\-\.]+)\s*step\s*=\s*([0-9]+)\s*"',
    re.IGNORECASE
)

def parse_title_time_step(dat_path: str):
    """Return (t_code: float, step: int) parsed from TITLE line."""
    with open(dat_path, "r", encoding="utf-8", errors="ignore") as f:
        head = "".join([next(f) for _ in range(5)])  # up to 5 lines
    m = TITLE_RE.search(head)
    if not m:
        raise ValueError(f"Could not parse TITLE line in {dat_path!r}. First lines were:\n{head}")
    t_code = float(m.group(1))
    step = int(m.group(2))
    return t_code, step

def code_time_to_physical(t_code: float,
                          start_iso: str = START_ISO) -> datetime:
    """Convert code-time 't_code' to physical UTC datetime using t_phys = t_code * (AU / V)."""
    start_dt = datetime.fromisoformat(start_iso)
    time_sec = t_code
    return start_dt + timedelta(seconds=time_sec)



def _parse_variables(line: str) -> List[str]:
    return re.findall(r'"([^"]+)"', line)

def _parse_zone_header(line: str) -> Dict[str, Any]:
    """Extract key=val pairs, preserving bracketed/quoted values."""
    props = {}
    raw = line.split('ZONE', 1)[1]
    tokens = re.findall(r'(\w+)=(".*?"|\([^)]*\)|\[.*?\]|[^,\s]+)', raw)
    for k, v in tokens:
        v = v.strip().strip('"')
        props[k] = v
    if 'I' in props: props['I'] = int(props['I'])
    if 'J' in props: props['J'] = int(props['J'])
    if 'N' in props: props['N'] = int(props['N'])
    if 'E' in props: props['E'] = int(props['E'])
    return props

def rotation_angle_radians(phys_dt: datetime,
                           rot_rate_deg_per_day: float = 14.184,
                           t0_iso: str = START_ISO) -> float:
    """Angle Sun has rotated since t0, in radians (CCW positive)."""
    t0 = datetime.fromisoformat(t0_iso)
    dt_days = (phys_dt - t0).total_seconds() / 86400.0
    return np.deg2rad(rot_rate_deg_per_day * dt_days)

def rotate_xy(x: np.ndarray, y: np.ndarray, theta: float) -> Tuple[np.ndarray, np.ndarray]:
    c, s = np.cos(theta), np.sin(theta)
    return c*x - s*y, s*x + c*y

def rotate_components_xy(cx: np.ndarray, cy: np.ndarray, theta: float) -> Tuple[np.ndarray, np.ndarray]:
    return rotate_xy(cx, cy, theta)

# ---------------------------
# Reader
# ---------------------------

def read_tecplot_slice_ascii(filename: str) -> Dict[str, Any]:
    """
    Reads Tecplot ASCII slice (BLOCK).
    Returns:
      {
        'variables': [..],
        'zones': [
           # Structured zone:
           {'kind':'STRUCT', 'I':int, 'J':int, 'props':dict, 'name':str,
            'blocks':{
                'node_blocks':[X_flat,Y_flat,Z_flat],
                'cell_blocks':{varname: flat array}
            }},
           # FE zone (ZONETYPE=FEQUADRILATERAL or N/E without I/J):
           {'kind':'FE', 'N':int, 'E':int, 'props':dict, 'name':str,
            'nodes':(xN,yN,zN),                 # each length N
            'cells':{varname: np.ndarray len E},# cell-centered vars
            'conn': np.ndarray shape (E,4),     # 1-based indices
           },
           ...
        ]
      }
    """
    import numpy as np, re

    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip()]

    variables: Optional[List[str]] = None
    zones: List[Dict[str, Any]] = []
    i = 0

    name_counts: Dict[str, int] = {}

    def unique_name(t_raw: str) -> str:
        c = name_counts.get(t_raw, 0); name_counts[t_raw] = c + 1
        return f"{t_raw}#{c}"

    def parse_numbers(s: str) -> List[float]:
        # grab numeric tokens (floats/ints/scientific); commas/spaces ok
        return [float(x) for x in re.findall(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?', s)]

    while i < len(lines):
        s = lines[i]
        su = s.upper().lstrip()

        if su.startswith('VARIABLES'):
            variables = _parse_variables(s)
            i += 1
            continue

        if su.startswith('TITLE') or su.startswith('DATASETAUXDATA'):
            i += 1
            continue

        if su.startswith('ZONE'):
            props = _parse_zone_header(s)
            t_raw = str(props.get('T','zone'))
            zonetype = props.get('ZONETYPE','').upper()
            has_IJ = ('I' in props and 'J' in props)
            has_NE = ('N' in props and 'E' in props)

            # -------- FE ZONE (deterministic read) --------
            if (not has_IJ and has_NE) or zonetype.startswith('FE'):
                n = int(props.get('N', 0)); e = int(props.get('E', 0))
                if variables is None:
                    raise ValueError("VARIABLES must appear before FE zones.")

                need_vals = 3*n + max(len(variables)-3, 0)*e  # floats to read first
                floats: List[float] = []
                ints: List[int] = []

                # consume following lines until next ZONE; first collect EXACTLY need_vals floats,
                # then collect connectivity ints (4*E expected)
                i += 1
                while i < len(lines) and not lines[i].upper().lstrip().startswith('ZONE'):
                    st = lines[i].lstrip()
                    if st.startswith('#'):
                        i += 1; continue
                    nums = re.findall(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?', lines[i])
                    if nums:
                        # how many floats still needed?
                        remain = need_vals - len(floats)
                        if remain > 0:
                            take = min(remain, len(nums))
                            # first 'take' tokens go into floats
                            floats.extend(float(x) for x in nums[:take])
                            # leftover tokens on the same line are connectivity ints
                            if take < len(nums):
                                ints.extend(int(float(x)) for x in nums[take:])
                        else:
                            # all tokens now belong to connectivity
                            ints.extend(int(float(x)) for x in nums)
                    i += 1

                if len(floats) != need_vals:
                    raise ValueError(f"{t_raw}: FE float block len {len(floats)} != expected {need_vals} (3*N + (V-3)*E).")

                # split float block into nodes + cell values
                arr = np.asarray(floats, dtype=float)
                xN = arr[0:n]; yN = arr[n:2*n]; zN = arr[2*n:3*n]
                cells = {}
                off = 3*n
                for vname in variables[3:]:
                    cells[vname] = arr[off:off+e]; off += e

                # connectivity
                if len(ints) < 4*e:
                    raise ValueError(f"{t_raw}: FE connectivity length {len(ints)} < 4*E={4*e}.")
                conn = np.asarray(ints[:4*e], dtype=int).reshape(e,4)

                zones.append({
                    'kind':'FE',
                    'N': n, 'E': e,
                    'props': props,
                    'name': unique_name(t_raw),
                    'nodes': (xN, yN, zN),
                    'cells': cells,
                    'conn': conn
                })
                # do not i+=1 here; outer loop continues with current i (already at next ZONE or EOF)
                continue

            # -------- STRUCTURED ZONE --------
            if not has_IJ:
                raise ValueError(f"ZONE missing I/J: {s}")

            I, J = int(props['I']), int(props['J'])
            if I < 2 or J < 2:
                # fast-forward to next ZONE
                i += 1
                while i < len(lines) and not lines[i].upper().lstrip().startswith('ZONE'):
                    i += 1
                continue

            varloc = props.get('VARLOCATION','')
            if varloc and ('CELLCENTERED' not in varloc.upper()):
                raise ValueError(f"Unsupported VARLOCATION for zone: {varloc}. Expected cell-centered vars 4..N.")

            # Read exactly: 3*(I*J) + (len(vars)-3)*((I-1)*(J-1)) floats for this zone
            n_nodes = I*J
            n_cells = (I-1)*(J-1) if variables and len(variables)>3 else 0
            need_vals = 3*n_nodes + max(len(variables)-3,0)*n_cells

            vals: List[float] = []
            i += 1
            while i < len(lines) and not lines[i].upper().lstrip().startswith('ZONE') and len(vals) < need_vals:
                st = lines[i].lstrip()
                if st.startswith('#'):
                    i += 1; continue
                vals.extend(parse_numbers(lines[i]))
                i += 1

            if len(vals) != need_vals:
                raise ValueError(
                    f"{t_raw}: STRUCT data len {len(vals)} != expected {need_vals} "
                    f"(3*{n_nodes} + {(len(variables)-3) if variables else 0}*{n_cells}) for I={I}, J={J}."
                )

            arr = np.asarray(vals, dtype=float)
            off = 0
            node_blocks = [arr[off:off+n_nodes]]; off += n_nodes
            node_blocks.append(arr[off:off+n_nodes]); off += n_nodes
            node_blocks.append(arr[off:off+n_nodes]); off += n_nodes

            cell_blocks_list = []
            for _ in range(len(variables)-3):
                cell_blocks_list.append(arr[off:off+n_cells]); off += n_cells

            zones.append({
                'kind':'STRUCT',
                'I': I, 'J': J,
                'props': props,
                'name': unique_name(t_raw),
                'blocks': {
                    'node_blocks': node_blocks,
                    'cell_blocks': dict(zip(variables[3:], cell_blocks_list))
                }
            })
            continue

        # skip stray non-ZONE lines here (the ZONE loops consume until next ZONE)
        i += 1

    if variables is None:
        raise ValueError("VARIABLES line not found.")

    return {'variables': variables, 'zones': zones}

# ---------------------------
# Helpers: AMR compositing & reshape
# ---------------------------

# ---------------------------
# Combined plot
# ---------------------------

def plot_with_geometry(
    geom: dict,
    varname: str,
    cell_order: str = "C",        # use "F" for p
    cmap: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    colorbar: bool = True,
    figsize=(10, 8),
    out_path: str = ".",
):
    """
    Render one variable using precomputed geometry (no recomputation of centers/triangulation).
    Requires: get_values_for_var(), UNIT_SCALE, RADIAL_SCALE, COLOR_SETTINGS, LABELS.
    Returns: path to saved PNG.
    """
    # Unpack geometry
    x_all = geom["x"]; y_all = geom["y"]; rc_all = geom["rc"]
    triang = geom["tri"]
    inner_radius = geom["inner"]; outer_radius = geom["outer"]
    step = geom["step"]; t_code = geom["t_code"]; phys_dt = geom["phys_dt"]; theta = geom["theta"]

    # Fetch cell values in the SAME order as geometry was assembled
    # if varname is T, calulate temperature from pressure and density
    if varname == "T":
        density = get_values_for_var(geom, "density", cell_order=cell_order, rotate_components=False)
        pressure = get_values_for_var(geom, "P", cell_order=cell_order, rotate_components=False)
        # Avoid division by zero
        v_flat = pressure/(2.0*(density/1.67262192e-24)*1.3806505e-16)
    else:
        v_flat = get_values_for_var(geom, varname, cell_order=cell_order, rotate_components=False)

    # Units → radial scaling
    v_raw = v_flat * UNIT_SCALE.get(varname, 1.0)
    pwr = RADIAL_SCALE.get(varname, 0.0)
    v_all = v_raw * ((rc_all/AU_CM) ** pwr) if pwr != 0 else v_raw

    # Colormap / range (allow overrides)
    cs = COLOR_SETTINGS.get(varname, {})
    cmap = cmap or cs.get("cmap", "plasma")
    vmin = cs.get("vmin", vmin)
    vmax = cs.get("vmax", vmax)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    if vmin is not None and vmax is not None:
        levels = np.linspace(vmin, vmax, 256)
        tpc = ax.tricontourf(triang, v_all, levels=levels, cmap=cmap, extend='both')
    else:
        tpc = ax.tricontourf(triang, v_all, levels=256, cmap=cmap, extend='both')

    if colorbar:
        cbar = fig.colorbar(tpc, ax=ax)
        cbar.set_label(LABELS.get(varname, varname))
        if vmin is not None and vmax is not None:
            ticks = np.linspace(vmin, vmax, 11)
            cbar.set_ticks(ticks)
            # scientific for tiny/huge ranges or density/temperature
            if varname in ["density", "T"] or (abs(vmax) < 1e-2 or abs(vmax) > 1e4):
                cbar.ax.set_yticklabels([f"{t:.1e}" for t in ticks])
            else:
                cbar.ax.set_yticklabels([f"{t:.2f}" for t in ticks])

    ax.xaxis.set_major_formatter(FuncFormatter(cm_to_AU))
    ax.yaxis.set_major_formatter(FuncFormatter(cm_to_AU))
    ax.set_xlabel("X (AU)"); ax.set_ylabel("Y (AU)")
    ax.set_aspect("equal", adjustable="box")
    # If you want hard limits:
    # ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2)

    # Title (step padded, 7 sig figs)
    title_left = f"{LABELS.get(varname, varname)}   R=[{inner_radius/AU_CM:.2f}, {outer_radius/AU_CM:.2f}] AU"
    try:
        step_str = f"{step:06d}" if (step is not None and step >= 0) else "------"
        t_str = f"{t_code:.7g}" if t_code is not None else "----"
        utc_str = phys_dt.strftime('%Y-%m-%d %H:%M:%S') if phys_dt else "----"
        rot_str = f"rot={np.rad2deg(-theta):.2f}°" if theta != 0 else "rot=0°"
        ax.set_title(f"{title_left}\nUTC: {utc_str}   {rot_str}")
    except Exception:
        ax.set_title(title_left)

    # Clean inner hole
    ax.add_patch(Circle((0, 0), inner_radius, facecolor="white", edgecolor="none", zorder=10))

    fig.tight_layout()
    out_filename = f"{out_path.rstrip('/')}/{varname}_{(step if step and step>=0 else 0):06d}.png"
    fig.savefig(out_filename, dpi=150)
    plt.close(fig)
    return out_filename

# ---------------------------
# Main: parallel plot slices
# ---------------------------

# =========================
# Parallel per-file (simple)
# =========================
import os, re
from multiprocessing import Pool, cpu_count
import matplotlib
matplotlib.use("Agg")  # safe for subprocesses with Matplotlib

VARIABLES = ["density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp", "T"]

def _render_one_file(args):
    filepath, out_dir = args
    # Reuse the functions already loaded in this script
    from __main__ import prepare_geometry, plot_with_geometry
    # Build geometry ONCE for this file
    geom = prepare_geometry(
        filename=filepath,
        reshape_order="C",
        corotate=COROTATE,                          # inertial? -> set False if you want
        rot_rate_deg_per_day=14.1843971631,
    )
    # Render all variables quickly (no geometry recompute)
    outputs = []
    for var in VARIABLES:
        out_png = plot_with_geometry(
            geom,
            varname=var,
            cell_order="C",
            out_path=out_dir,
        )
        outputs.append((os.path.basename(filepath), var, os.path.basename(out_png)))
    return outputs

if __name__ == "__main__":
    
    os.makedirs(out_dir, exist_ok=True)

    MIN_ITER = -1

    #Search all density*.png files in out_dir to determine max iteration already done
    if not PROCESS_ALL:
        existing_pngs = [f for f in os.listdir(out_dir) if re.match(r'density_\d+\.png$', f)]
        max_done_iter = -1
        for f in existing_pngs:
            m = re.search(r'density_(\d+)\.png$', f)
            if m:
                iter_num = int(m.group(1))
                if iter_num > max_done_iter:
                    max_done_iter = iter_num
        MIN_ITER = max_done_iter
        print(f"[INFO] Skipping files with iter <= {MIN_ITER} based on existing density PNGs.")

    all_files = [f for f in os.listdir(slice_folder) if re.match(r'.*\.z\.\d+\.dat$', f)]
    if not all_files:
        raise FileNotFoundError(f"No slice files found in {slice_folder!r} matching '*.z.*.dat'")

    def iter_num(fname: str) -> int:
        return int(re.search(r'\.z\.(\d+)\.dat$', fname).group(1))

    slice_files = sorted((f for f in all_files if iter_num(f) > MIN_ITER), key=iter_num)
    if not slice_files:
        raise SystemExit(f"No files with iter > {MIN_ITER} found in {slice_folder}")

    tasks = [(os.path.join(slice_folder, f), out_dir) for f in slice_files]

    # Use N-1 workers by default; change to a fixed number (e.g., 8) if you prefer
    # workers = max(1, (cpu_count() or 2) - 1)
    workers = 8
    print(f"[INFO] Processing {len(slice_files)} files × {len(VARIABLES)} vars using {workers} workers...")

    done = 0
    total = len(slice_files) * len(VARIABLES)
    with Pool(processes=workers) as pool:
        for file_results in pool.imap_unordered(_render_one_file, tasks, chunksize=1):
            # file_results is a list of (file, var, out_png) for that file
            done += len(file_results)
            file_shown = file_results[0][0] if file_results else "?"
            print(f"[{done}/{total}] {file_shown}: {len(file_results)} plots")

    print("[INFO] All plots complete.")
            