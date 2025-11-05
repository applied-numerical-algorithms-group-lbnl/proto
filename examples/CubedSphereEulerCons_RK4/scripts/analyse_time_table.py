import re
import pandas as pd
import plotly.express as px

# ------- CONFIG -------

FILENAME = "/Users/talwindersingh/Desktop/My_Computer/Work/UAH/Current_projects/SWQU/proto/examples/CubedSphereEulerCons_RK4/exec_Convergence/45_DIM3_30_45_CubeSphereTest.time.table"
THRESHOLD_PCT = 0.00      # plot nodes with >= this % of total time.
VIEW_MODE = "treemap"   # "sunburst" or "treemap" #"sunburst" doesn't work well sometimes.

# ------- PARSER -------

def parse_gprof_callgraph(text: str):
    """
    Parse a gprof-like callgraph report into:
      - nodes: dict[id] = {id, name, time, calls}
      - edges: list of (parent_id, child_id)

    This version:
      * Skips everything before the FIRST
            ---------------------------------------------------------
        and starts parsing from the line below (which in your file
        is `[0]root ...`).
    """

    # 1) Skip everything before the first long line of dashes
    sep = "---------------------------------------------------------"
    pos = text.find(sep)
    if pos == -1:
        raise RuntimeError("Could not find the '---------------------------------------------------------' separator in file.")

    # Start parsing from just after that separator
    text = text[pos + len(sep):]

    # 2) Now split the remaining text into blocks by the same separator
    blocks = text.strip().split(sep)
    nodes = {}
    edges = []

    # Header line pattern for callgraph blocks:
    #   [0]root 316.62062 1
    #   [1] MMBEuler 315.12823 1
    header_re = re.compile(r'\[(\d+)\]\s*(.+?)\s+([0-9.]+)\s+(\d+)\s*$')

    for block in blocks:
        block = block.strip()
        if not block:
            continue

        lines = block.splitlines()
        if not lines:
            continue

        # Header line must look like: [id] name time calls
        m = header_re.match(lines[0])
        if not m:
            continue

        node_id = int(m.group(1))
        name = m.group(2).strip()
        total_time = float(m.group(3))
        calls = int(m.group(4))

        nodes[node_id] = {
            "id": node_id,
            "name": name,
            "time": total_time,
            "calls": calls,
        }

        # Child lines: "  94.5% 297.7679   10 ChildName [2]"
        for line in lines[1:]:
            line = line.rstrip()
            if "Total" in line:
                continue
            cm = re.match(r'\s*([\d.]+)%\s+([0-9.]+)\s+(\d+)\s+(.+?)\s+\[(\d+)\]', line)
            if not cm:
                continue

            parent_id = node_id
            child_id = int(cm.group(5))
            edges.append((parent_id, child_id))

            # Ensure child exists in nodes dict, even if its own header appears later
            if child_id not in nodes:
                nodes[child_id] = {
                    "id": child_id,
                    "name": cm.group(4).strip(),
                    "time": None,
                    "calls": int(cm.group(3)),
                }

    return nodes, edges


# ------- BUILD DF + FILTERING -------

def build_sunburst_df(nodes, edges, threshold_pct=0.1):
    """
    Build DataFrame for Plotly (sunburst or treemap).

    - Keeps nodes whose global % of total time >= threshold_pct
      (based on time in their header line),
      BUT also keeps all ancestors of any kept node, so hierarchy is intact.
    """
    # total time from root [0]
    if 0 in nodes and nodes[0]["time"] is not None:
        total_time = nodes[0]["time"]
    else:
        total_time = max(nd["time"] for nd in nodes.values() if nd["time"] is not None)

    # parent map: first parent per child
    parent = {}
    for p, c in edges:
        if c not in parent:
            parent[c] = p
    parent[0] = ""  # root has no parent

    # compute pct_total per node
    pct_by_id = {}
    for nid, nd in nodes.items():
        value = nd["time"] or 0.0
        pct = 100.0 * value / total_time if total_time > 0 else 0.0
        pct_by_id[nid] = pct

    # 1) initial keep set by threshold
    keep_ids = set()
    for nid, pct in pct_by_id.items():
        if nid == 0 or pct >= threshold_pct:
            keep_ids.add(nid)

    # 2) closure under ancestry: include all parents of kept nodes
    added = True
    while added:
        added = False
        for nid in list(keep_ids):
            p = parent.get(nid, "")
            if isinstance(p, int) and p not in keep_ids:
                keep_ids.add(p)
                added = True

    # Build rows for DataFrame
    rows = []
    for nid, nd in nodes.items():
        if nid not in keep_ids:
            continue
        value = nd["time"] or 0.0
        pct_total = pct_by_id[nid]
        rows.append({
            "id": nid,
            "label": f"{nd['name']} ({pct_total:.2f}%)",
            "parent": parent.get(nid, ""),
            "value": value,
            "pct_total": pct_total,
        })

    return pd.DataFrame(rows)


# ------- MAIN -------

def main():
    with open(FILENAME, "r") as f:
        text = f.read()

    nodes, edges = parse_gprof_callgraph(text)
    df = build_sunburst_df(nodes, edges, threshold_pct=THRESHOLD_PCT)

    print(f"Total nodes parsed: {len(nodes)}")
    print(f"Nodes shown in plot: {len(df)} (threshold = {THRESHOLD_PCT}%)")
    print(df.head())

    if VIEW_MODE == "sunburst":
        fig = px.sunburst(
            df,
            ids="id",
            names="label",
            parents="parent",
            values="value",
            color="pct_total",
            color_continuous_scale="Viridis",
            branchvalues="total",
        )
    elif VIEW_MODE == "treemap":
        fig = px.treemap(
            df,
            ids="id",
            names="label",
            parents="parent",
            values="value",
            color="pct_total",
            color_continuous_scale="Viridis",
        )
    else:
        raise ValueError(f"Unknown VIEW_MODE = {VIEW_MODE}")

    # 🔹 Smaller, cleaner hover box
    hover_tmpl = (
        "<b>%{label}</b><br>"           # label already has function name + %
        "time = %{value:.3f} s<br>"
        "% of total = %{color:.2f}%<br>"
        "id = %{id}<br>"
        "parent = %{parent}<br>"
        "<extra></extra>"               # removes the grey trace name at bottom
    )
    fig.update_traces(hovertemplate=hover_tmpl)

    # 🔹 Smaller hover font so it doesn’t feel huge when zoomed/focused
    fig.update_layout(
        hoverlabel=dict(font_size=10),
        title=f"Callgraph Time Breakdown (≥ {THRESHOLD_PCT}% total, view={VIEW_MODE})",
        margin=dict(t=40, l=0, r=0, b=0),
    )

    fig.show()


if __name__ == "__main__":
    main()