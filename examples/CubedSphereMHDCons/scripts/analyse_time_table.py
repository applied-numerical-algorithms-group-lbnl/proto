import re
import pandas as pd
import plotly.express as px
import numpy as np
import colorsys

# ---------------- CONFIG ----------------

FILENAME = "/Users/talwindersingh/Desktop/My_Computer/Work/UAH/Current_projects/SWQU/proto/examples/CubedSphereMHDCons/exec/180_DIM3_15_45_CubeSphereTest.time.table"
# FILENAME = "/Users/talwindersingh/Desktop/My_Computer/Work/UAH/Current_projects/SWQU/proto/examples/CubedSphereMHDCons/exec/90_DIM3_30_90_CubeSphereTest.time.table"
THRESHOLD_PCT = 0.1   # keep nodes with >= this % of total time
VIEW_MODE = "treemap"  # only treemap for now


# ------------- STEP 1: EXTRACT TIMER REPORT BLOCK ------------- #

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

    # 4) Now collect timer lines [NN] ... until we hit the long dashed
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
            # This is a timer line
            timer_lines.append((idx + 1, line))  # store 1-based line number
            seen_any_timer = True

        idx += 1

    return timer_lines


# ------------- STEP 2: COMPUTE LINEAGE BY INDENTATION ------------- #

def compute_lineage_for_lines(lines_with_numbers):
    """
    Given a list of (line_number, line_text) for the timer tree,
    compute lineage of each line based on indentation.

    Returns: dict[line_number] = [line_number, parent, grandparent, ...]
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


# ------------- STEP 3: PARSE TIMER FIELDS FROM EACH LINE ------------- #

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
    # Example line:
    # [0] root 493.88826 100.0% 1 8954189530   18.1 MFlops
    timer_re = re.compile(
        r'^\s*\[(\d+)\]\s+(.+?)\s+([0-9.]+)\s+([0-9.]+)%\s+(\d+)\b'
    )

    records = {}
    for lnum, line in lines_with_numbers:
        m = timer_re.match(line)
        if not m:
            # Shouldn't happen if extract_timer_lines_with_numbers worked,
            # but we'll just skip if it does.
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


# ------------- STEP 4: BUILD DATAFRAME FOR TREEMAP ------------- #

def build_treemap_df(records, lineage, parent_map, threshold_pct=0.0):
    """
    records: dict[line] = parsed timer info
    lineage: dict[line] = [self, parent, grandparent, ...]
    parent_map: dict[line] = parent_line or None

    Returns a DataFrame with columns:
      id, label, parent, value, pct_total, depth, timer_id, name
    suitable for Plotly treemap.
    """

    # Identify roots: lines that have parent None
    roots = [ln for ln, p in parent_map.items() if p is None]

    # Use the first root for total-time reference (usually [0] root)
    if roots:
        root_line = roots[0]
    else:
        # fallback: pick smallest line number
        root_line = min(records.keys())

    # total time from root (just for sanity; pct is already absolute)
    total_time = records[root_line]["time"]

    # initial keep set: nodes with pct >= threshold
    keep_lines = set()
    for lnum, rec in records.items():
        if rec["pct"] >= threshold_pct or lnum == root_line:
            keep_lines.add(lnum)

    # closure under ancestry: include all ancestors of kept lines
    added = True
    while added:
        added = False
        for lnum in list(keep_lines):
            chain = lineage.get(lnum, [])
            # chain[0] is self; ancestors are chain[1:]
            for anc in chain[1:]:
                if anc not in keep_lines:
                    keep_lines.add(anc)
                    added = True

    # depth = length of lineage - 1
    def get_depth(lnum):
        chain = lineage.get(lnum, [lnum])
        return len(chain) - 1

    rows = []
    for lnum, rec in records.items():
        if lnum not in keep_lines:
            continue

        pct_total = rec["pct"]
        value = rec["time"]
        depth = get_depth(lnum)

        parent_line = parent_map.get(lnum)
        parent_id = parent_line if parent_line is not None else ""

        label = f"{rec['name']} ({pct_total:.2f}%)"

        rows.append({
            "id": lnum,               # use file line number as node id
            "label": label,
            "parent": parent_id,      # parent is also a line number (or "")
            "value": value,           # time in seconds
            "pct_total": pct_total,
            "depth": depth,
            "timer_id": rec["timer_id"],
            "name": rec["name"],
            "calls": rec["calls"],
        })

    return pd.DataFrame(rows)


# ------------- STEP 5: MAIN + PLOTTING ------------- #

def main():


    # Match timer lines and expected "strict" format (id, name, time, %, calls)
    timer_line_re = re.compile(r'^\s*\[(\d+)\]')
    strict_re = re.compile(
        r'^\s*\[(\d+)\]\s+(.+?)\s+([0-9.]+)\s+([0-9.]+)%\s+(\d+)\b'
    )

    malformed_lines = []
    with open(FILENAME, "r") as f:
        for lineno, line in enumerate(f, start=1):
            # Stop checking once the dashed separator appears
            if "------------------------------" in line:
                break
            
            
            # Check only lines that look like timers
            if timer_line_re.match(line) and not strict_re.match(line):
                print(f"⚠️ Malformed timer line at {lineno}: {line.rstrip()}")
                malformed_lines.append((lineno, line.rstrip()))
    if malformed_lines:
        print(f"\nTotal malformed timer lines found: {len(malformed_lines)}")
        return
                

    with open(FILENAME, "r") as f:
        text = f.read()

    # Extract timer tree lines
    timer_lines = extract_timer_lines_with_numbers(text)
    print(f"Found {len(timer_lines)} timer lines in Timer report 0 block.")

    # Compute lineage and parent map
    lineage, parent_map = compute_lineage_for_lines(timer_lines)

    # Parse fields (id, name, time, pct, calls)
    records = parse_timer_fields(timer_lines)

    # Build DataFrame for treemap
    df = build_treemap_df(records, lineage, parent_map, threshold_pct=THRESHOLD_PCT)

    # Save for inspection (if desired)
    # out_csv = FILENAME + f".timer_tree.pct{THRESHOLD_PCT:.2f}.csv"
    # df.to_csv(out_csv, index=False)
    # print(f"DataFrame saved to: {out_csv}")
    print(f"Nodes shown in plot: {len(df)}  (threshold = {THRESHOLD_PCT}%)")
    print(df.head())

    # --- RANDOM COLORS, BRIGHTER FOR DEEPER NODES ---
    rng = np.random.default_rng(123)
    max_depth = max(df["depth"].max(), 1)

    def depth_color(depth):
        h = rng.random()
        s = 0.6
        v = 0.4 + 0.6 * (depth / max_depth)
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        return f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})"

    df["color"] = df["depth"].apply(depth_color)

    if VIEW_MODE == "treemap":
        fig = px.treemap(
            df,
            ids="id",
            names="label",
            parents="parent",
            values="value",
            custom_data=["pct_total", "id", "parent", "timer_id", "name", "calls"],
        )
        fig.update_traces(marker=dict(colors=df["color"]))
    else:
        raise ValueError(f"Unknown VIEW_MODE = {VIEW_MODE}")

    # Hover template
    hover_tmpl = (
        "<b>%{label}</b><br>"            # name
        "parent: %{parent}<br>"
        "time = %{value:.6f} s<br>"
        "% of total = %{customdata[0]:.2f}%<br>"
        # "file line = %{customdata[1]}<br>"
        # "parent line = %{customdata[2]}<br>"
        # "timer_id = %{customdata[3]}<br>"
        "calls = %{customdata[5]}<br>"
        "<extra></extra>"
    )
    fig.update_traces(hovertemplate=hover_tmpl)
    file_without_path = FILENAME.split("/")[-1]
    fig.update_layout(
        hoverlabel=dict(font_size=10),
        title=f"Timer Tree Time Breakdown (≥ {THRESHOLD_PCT}% total, view={VIEW_MODE}) | File: {file_without_path}",
        margin=dict(t=40, l=0, r=0, b=0),
    )

    #save to html (if desired)
    # out_html = FILENAME + f".timer_tree.pct{THRESHOLD_PCT:.2f}.html"
    # fig.write_html(out_html)
    # print(f"Treemap plot saved to: {out_html}")

    fig.show()


if __name__ == "__main__":
    main()