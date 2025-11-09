import re
import pandas as pd
import plotly.express as px
import numpy as np
import colorsys

# ---------------- CONFIG ----------------

FILES = [
    "/Users/talwindersingh/Desktop/My_Computer/Work/UAH/Current_projects/SWQU/proto/examples/CubedSphereMHDCons/exec/180_DIM3_15_45_CubeSphereTest.time.table",
    "/Users/talwindersingh/Desktop/My_Computer/Work/UAH/Current_projects/SWQU/proto/examples/CubedSphereMHDCons/exec/180_DIM3_30_20_CubeSphereTest.time.table",
]
NGHOST = 5  # number of ghost cells in each direction

#Get workloads from filenames
def get_workload_from_filename(filename):
    m = re.search(r'DIM3_(\d+)_(\d+)_CubeSphereTest', filename)
    if not m:
        raise ValueError(f"Could not extract workload from filename: {filename}")
    n1 = int(m.group(1))
    n2 = int(m.group(2))
    return (n1 + NGHOST) * (n1 + NGHOST) * (n2 + NGHOST)

WORKLOADS = [get_workload_from_filename(f) for f in FILES]

THRESHOLD_PCT_REF_RUN = 0.1   # threshold in % of total (reference run)
REF_RUN_INDEX = 1             # which run defines treemap areas (0 or 1)
VIEW_MODE = "treemap"


# ---------- STEP 1: COMMON PARSING UTILITIES (same logic as working version) ---------- #

def extract_timer_lines_with_numbers(text: str):
    """
    From the full file text, extract the timer tree block under

        Timer report 0 (.... timers)
        --------------

    and before the first long dashed separator

        ---------------------------------------------------------

    Returns: list of (line_number, line_text) (1-based file line numbers)
    Only lines that look like timers (start with [NN]) are kept.
    """
    lines = text.splitlines()

    # 1) Find "Timer report 0" line
    start_idx = None
    for i, line in enumerate(lines):
        if "Timer report 0" in line:
            start_idx = i
            break

    if start_idx is None:
        raise RuntimeError("Could not find 'Timer report 0' in file.")

    # 2) Move to the line after "Timer report 0 (...)"
    idx = start_idx + 1

    # 3) Skip the short dashed separator(s) like "--------------"
    while idx < len(lines):
        s = lines[idx].strip()
        if s and set(s) == {"-"}:
            idx += 1
        else:
            break

    # 4) Collect timer lines [NN] ... until we hit the long dashed
    #    separator that starts the callgraph part
    timer_lines = []
    seen_any_timer = False
    timer_re = re.compile(r'^\s*\[\d+\]')  # line starts with [NN]

    while idx < len(lines):
        line = lines[idx]
        stripped = line.strip()

        # long dashed separator (20+ dashes) -> stop once we've seen timers
        if stripped and set(stripped) == {"-"} and len(stripped) >= 20:
            if seen_any_timer:
                break

        if timer_re.match(line):
            timer_lines.append((idx + 1, line))  # store 1-based line number
            seen_any_timer = True

        idx += 1

    return timer_lines


def compute_lineage_for_lines(lines_with_numbers):
    """
    Given a list of (line_number, line_text) for the timer tree,
    compute lineage of each line based on indentation.

    Returns:
      lineage: dict[line_number] = [line_number, parent, grandparent, ...]
      parent: dict[line_number]  = direct parent line_number or None
    """

    stack = []      # (indent_len, line_number)
    parent = {}     # parent[line_number] = parent_line_number or None

    for lnum, raw_line in lines_with_numbers:
        if raw_line.strip() == "":
            continue

        # Normalize tabs -> spaces, then count leading spaces as indent
        s = raw_line.expandtabs(4)
        indent_len = len(s) - len(s.lstrip(" "))

        # Pop until the top has strictly smaller indent
        while stack and stack[-1][0] >= indent_len:
            stack.pop()

        # Parent is now the last element in stack (if any)
        if stack:
            parent[lnum] = stack[-1][1]
        else:
            parent[lnum] = None  # top-level within timer tree

        # Push this line on the stack
        stack.append((indent_len, lnum))

    # Build lineage: for each line, walk up parent[]
    lineage = {}
    for lnum, raw_line in lines_with_numbers:
        if raw_line.strip() == "":
            continue
        chain = []
        cur = lnum
        while cur is not None:
            chain.append(cur)
            cur = parent.get(cur)
        lineage[lnum] = chain  # [self, parent, grandparent, ...]

    return lineage, parent


def parse_timer_fields(lines_with_numbers):
    """
    Parse each timer line:

        [id] name time pct% calls ...

    Returns: dict[line_number] = {
        "line": lnum,
        "timer_id": int,
        "name": str,
        "time": float (seconds),
        "pct": float,
        "calls": int,
    }
    """
    timer_re = re.compile(
        r'^\s*\[(\d+)\]\s+(.+?)\s+([0-9.]+)\s+([0-9.]+)%\s+(\d+)\b'
    )

    records = {}
    for lnum, line in lines_with_numbers:
        m = timer_re.match(line)
        if not m:
            continue

        timer_id = int(m.group(1))
        name     = m.group(2).strip()
        time_val = float(m.group(3))
        pct_val  = float(m.group(4))
        calls    = int(m.group(5))

        records[lnum] = {
            "line": lnum,
            "timer_id": timer_id,
            "name": name,
            "time": time_val,
            "pct": pct_val,
            "calls": calls,
        }

    return records


# ---------- STEP 2: BUILD PATH KEYS (function + full lineage of names) ---------- #

def build_path_dict(records, lineage):
    """
    records: dict[line] = {name, time, pct, calls, ...}
    lineage: dict[line] = [self_line, parent_line, grandparent_line, ...]

    Returns:
      paths: dict[path_tuple] = {
          "time": total_time_for_this_path,
          "pct": total_pct_for_this_path,
          "calls": total_calls_for_this_path,
          "sample_line": one_line_number_for_reference,
          "name": leaf_name,
      }

    where path_tuple = (root_name, ..., parent_name, name)
    """
    # map line -> name for quick lookup
    name_by_line = {ln: rec["name"] for ln, rec in records.items()}

    paths = {}

    for lnum, rec in records.items():
        chain = lineage.get(lnum, [lnum])   # [self, parent, ...]
        # build path from root -> leaf (reverse the order)
        name_chain = [name_by_line[i] for i in chain]
        path = tuple(reversed(name_chain))

        if path not in paths:
            paths[path] = {
                "time": 0.0,
                "pct": 0.0,
                "calls": 0,
                "sample_line": lnum,
                "name": name_chain[0],  # leaf name is chain[0] before reversing
            }

        paths[path]["time"]  += rec["time"]
        paths[path]["pct"]   += rec["pct"]
        paths[path]["calls"] += rec["calls"]

    return paths


# ---------- STEP 3: LOAD BOTH RUNS AND MATCH PATHS ---------- #

def load_run(filename):
    with open(filename, "r") as f:
        text = f.read()

    timer_lines = extract_timer_lines_with_numbers(text)
    lineage, parent_map = compute_lineage_for_lines(timer_lines)
    records = parse_timer_fields(timer_lines)

    # root total time: the line whose name is "root", or the one with max pct
    root_candidates = [rec for rec in records.values() if rec["name"] == "root"]
    if root_candidates:
        root_time = root_candidates[0]["time"]
    else:
        root_time = max(rec["time"] for rec in records.values())

    paths = build_path_dict(records, lineage)

    return {
        "filename": filename,
        "timer_lines": timer_lines,
        "lineage": lineage,
        "parent_map": parent_map,
        "records": records,
        "paths": paths,
        "root_time": root_time,
    }


def combine_two_runs(run1, run2, W1, W2, threshold_pct_ref=0.0, ref_run_index=1):
    """
    run1, run2: results from load_run()
    W1, W2: workloads for run1 and run2
    threshold_pct_ref: threshold on pct_total of reference run
    ref_run_index: which run is reference for treemap (0 or 1)
    """

    # choose reference run
    runs = [run1, run2]
    Ws   = [W1, W2]
    ref  = runs[ref_run_index]
    ref_W = Ws[ref_run_index]
    other_idx = 1 - ref_run_index
    other = runs[other_idx]
    other_W = Ws[other_idx]

    # Paths dicts
    paths1 = run1["paths"]
    paths2 = run2["paths"]

    # intersection of paths
    common_paths = sorted(set(paths1.keys()) & set(paths2.keys()))

    rows = []
    for path in common_paths:
        p1 = paths1[path]
        p2 = paths2[path]

        t1 = p1["time"]
        t2 = p2["time"]

        # percent-of-total in each run
        pct1 = 100.0 * t1 / run1["root_time"] if run1["root_time"] > 0 else 0.0
        pct2 = 100.0 * t2 / run2["root_time"] if run2["root_time"] > 0 else 0.0

        # apply threshold on reference run's pct
        ref_pct = pct2 if ref_run_index == 1 else pct1
        if ref_pct < threshold_pct_ref:
            continue

        # scaling exponent p, if W1 != W2 and times > 0
        # p = None
        # if W1 != W2 and t1 > 0 and t2 > 0:
        #     ratio_t = t2 / t1
        #     ratio_W = W2 / W1
        #     if ratio_W != 1.0:
        #         p = np.log(ratio_t) / np.log(ratio_W)
        p = W1/W2 * (t2/t1)
        # label the exponent nicely
        p_display = p if p is not None else np.nan

        # create IDs from path (string)
        path_str = " > ".join(path)
        parent_path = path[:-1]
        parent_id = " > ".join(parent_path) if parent_path else ""

        leaf_name = path[-1]

        rows.append({
            "id": path_str,
            "parent": parent_id,
            "name": leaf_name,
            "full_path": path_str,
            "t1": t1,
            "t2": t2,
            "pct1": pct1,
            "pct2": pct2,
            "W1": W1,
            "W2": W2,
            "p": p_display,
        })

    # Ensure all parents exist as nodes (closure under ancestry)
    # Parents that do not exist in rows are created as zero-time "internal" nodes.
    existing_ids = {r["id"] for r in rows}
    parent_ids_needed = set()
    for r in rows:
        if r["parent"]:
            parent_ids_needed.add(r["parent"])

    # iterative closure
    added = True
    while added:
        added = False
        new_parents = set()
        for pid in list(parent_ids_needed):
            if pid not in existing_ids:
                # create synthetic parent node
                parent_path = pid.split(" > ")
                parent_parent = " > ".join(parent_path[:-1]) if len(parent_path) > 1 else ""
                leaf_name = parent_path[-1]
                rows.append({
                    "id": pid,
                    "parent": parent_parent,
                    "name": leaf_name,
                    "full_path": pid,
                    "t1": 0.0,
                    "t2": 0.0,
                    "pct1": 0.0,
                    "pct2": 0.0,
                    "W1": W1,
                    "W2": W2,
                    "p": np.nan,
                })
                existing_ids.add(pid)
                if parent_parent:
                    new_parents.add(parent_parent)
                added = True
        parent_ids_needed = new_parents

    return pd.DataFrame(rows)


# ---------- STEP 4: MAIN + PLOT ---------- #

def main():
    if len(FILES) != 2 or len(WORKLOADS) != 2:
        raise ValueError("This script currently expects exactly 2 FILES and 2 WORKLOADS.")

    run_data = [load_run(FILES[0]), load_run(FILES[1])]

    df = combine_two_runs(run_data[0], run_data[1],
                          WORKLOADS[0], WORKLOADS[1],
                          threshold_pct_ref=THRESHOLD_PCT_REF_RUN,
                          ref_run_index=REF_RUN_INDEX)

    print(f"Matched call-paths in both runs (after threshold): {len(df)}")
    # Save DataFrame for inspection (if desired)
    # out_csv = "scaling_treemap_data.csv"
    # df.to_csv(out_csv, index=False)
    # print(f"Scaling DataFrame saved to: {out_csv}")
    print(df.head())

    # depth from number of " > " segments
    df["depth"] = df["full_path"].apply(lambda s: s.count(">"))
    df["p_clamped"] = df["p"].clip(lower=0.0, upper=1.0)

    # # --- RANDOM COLORS BY DEPTH ---
    # rng = np.random.default_rng(123)
    # max_depth = max(df["depth"].max(), 1)

    # def depth_color(depth):
    #     h = rng.random()
    #     s = 0.6
    #     v = 0.4 + 0.6 * (depth / max_depth)
    #     r, g, b = colorsys.hsv_to_rgb(h, s, v)
    #     return f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})"

    # df["color"] = df["depth"].apply(depth_color)

    # Choose which run defines area: REF_RUN_INDEX
    value_col = "t2" if REF_RUN_INDEX == 1 else "t1"

    if VIEW_MODE == "treemap":
        # fig = px.treemap(
        #     df,
        #     ids="id",
        #     parents="parent",
        #     names="name",
        #     values=value_col,
        #     custom_data=["full_path", "t1", "t2", "pct1", "pct2", "W1", "W2", "p"],
        # )
        fig = px.treemap(
            df,
            ids="id",
            parents="parent",
            names="name",
            values=value_col,
            color="p_clamped",
            color_continuous_scale=[(0.0, "red"), (1.0, "green")],
            range_color=(0.0, 1.0),
            custom_data=["full_path", "t1", "t2", "pct1", "pct2", "W1", "W2", "p"],
        )
        fig.update_coloraxes(colorbar_title_text="p")
        # fig.update_traces(marker=dict(colors=df["color"]))
    else:
        raise ValueError(f"Unknown VIEW_MODE = {VIEW_MODE}")

    hover_tmpl = (
        # "<b>%{customdata[0]}</b><br>"           # full_path
        "name: %{label}<br>"
        "parent: %{parent}<br>"
        "t1 = %{customdata[1]:.6f} s<br>"
        "t2 = %{customdata[2]:.6f} s<br>"
        "pct1 = %{customdata[3]:.2f}%<br>"
        "pct2 = %{customdata[4]:.2f}%<br>"
        "W1 = %{customdata[5]}<br>"
        "W2 = %{customdata[6]}<br>"
        "p ((W1/W2) * (t2/t1)) = %{customdata[7]:.3f}<br>"
        "<extra></extra>"
    )
    fig.update_traces(hovertemplate=hover_tmpl)
    file1_without_path = FILES[0].split("/")[-1]
    file2_without_path = FILES[1].split("/")[-1]
    fig.update_layout(
        hoverlabel=dict(font_size=10),
        title=(
            f"Per-call-path Scaling Treemap | W1={WORKLOADS[0]}, W2={WORKLOADS[1]} | areas = t{REF_RUN_INDEX+1}<br>"
            f"File1: {file1_without_path} | File2: {file2_without_path}<br>"
        ),
        margin=dict(t=60, l=0, r=0, b=0),
    )

    #Save to html (if desired)
    out_html = "scaling_treemap.html"
    fig.write_html(out_html)
    print(f"Scaling treemap plot saved to: {out_html}")

    fig.show()


if __name__ == "__main__":
    main()