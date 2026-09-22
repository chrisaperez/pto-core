"""Graphic-spec model and exporters for the pto-core diagram suite (stdlib only).

One typed model (pages -> nodes, edges) is the single source of truth. It is exported to:
  * diagram-spec.json          machine-readable spec
  * pto-core-diagrams.drawio   multi-page draw.io / diagrams.net file (containers are real parents)
  * svg/<page>.svg             per-page SVG (upload to Canva, or open in any browser)
  * graphic-spec.md            human-readable node and edge tables

and checked for: unique ids, page bounds, containment, sibling overlap, dangling edges, an
estimate of whether each label fits its box, and WCAG contrast of every palette pair.
"""
from __future__ import annotations

import html
import json
from dataclasses import asdict, dataclass, field

PAGE_W, PAGE_H = 1920, 1080
FONT_SANS = "Inter, Helvetica, Arial, sans-serif"
FONT_MONO = "JetBrains Mono, Menlo, Consolas, monospace"

# ------------------------------------------------------------------ palette
KINDS = {
    #  kind        fill       stroke     text       dash   drawio shape hint            meaning
    "cpp":      ("#E3F0FD", "#1565C0", "#0D3B72", False, "rounded",  "C++20 engine or core component"),
    "py":       ("#E6F4EA", "#2E7D32", "#1B4D20", False, "rounded",  "Python orchestration (stdlib only)"),
    "build":    ("#EFE6FA", "#6A3FA0", "#3F2266", False, "rounded",  "build, packaging, CI, deploy"),
    "data":     ("#FFF4DB", "#B7791F", "#6B4A0E", False, "document", "artifact, format, stream"),
    "sync":     ("#DDF5F1", "#00796B", "#004D40", False, "process",  "concurrency primitive: ring, pool, lock"),
    "gate":     ("#FDE8E8", "#C62828", "#7F1D1D", False, "rounded",  "security gate"),
    "ext":      ("#F4F4F4", "#757575", "#424242", True,  "rounded",  "outside this repository (contract only)"),
    "panel":    ("#FFFFFF", "#B0BEC5", "#37474F", False, "container", "neutral grouping container"),
    "zone_cpp": ("#F5F9FF", "#1565C0", "#0D3B72", False, "container", "container: C++ engines"),
    "zone_py":  ("#F4FAF5", "#2E7D32", "#1B4D20", False, "container", "container: Python layer"),
    "zone_gate": ("#FFF6F6", "#C62828", "#7F1D1D", True, "container", "container: trust boundary"),
    "boundary": ("#EFE6FA", "#6A3FA0", "#3F2266", False, "band",      "the process boundary"),
    "callout":  ("#FFFBEA", "#E0B84C", "#5B4300", False, "rounded",  "callout: a finding worth reading"),
    "note":     ("#FFFFFF", "#FFFFFF", "#37474F", False, "text",      "annotation, no box"),
}
CONTAINER_KINDS = {"panel", "zone_cpp", "zone_py", "zone_gate"}

EDGES = {
    #  kind        colour     width dash   meaning
    "payload":  ("#37474F", 2.5, None,  "bytes, buffer handles, artifacts"),
    "control":  ("#1565C0", 1.6, "6,4", "call, argv, configuration"),
    "error":    ("#C62828", 1.6, "2,4", "abort, refusal, failure exit"),
    "contract": ("#6A3FA0", 3.0, None,  "crosses the process boundary"),
}


# ------------------------------------------------------------------ model
@dataclass
class Node:
    id: str
    kind: str
    x: int
    y: int
    w: int
    h: int
    label: str = ""
    body: list = field(default_factory=list)
    parent: str | None = None
    mono: bool = False
    table: dict | None = None      # {"cols": [px..], "header": [..], "rows": [[cell_lines..]..]}
    size: int = 13                 # body font size


@dataclass
class Edge:
    id: str
    src: str
    dst: str
    kind: str = "payload"
    label: str = ""
    s: str = "r"                   # exit side of src:  l r t b
    d: str = "l"                   # entry side of dst
    via: list = field(default_factory=list)
    s_at: float = 0.5              # where along the exit side (0..1)
    d_at: float = 0.5              # where along the entry side (0..1)


@dataclass
class Page:
    id: str
    name: str
    title: str
    subtitle: str
    footer: str
    nodes: list = field(default_factory=list)
    edges: list = field(default_factory=list)
    talk: list = field(default_factory=list)


class PageBuilder:
    def __init__(self, pid, name, title, subtitle, footer):
        self.page = Page(pid, name, title, subtitle, footer)
        self._n = 0

    def box(self, id, kind, x, y, w, h, label="", body=(), parent=None, **kw):
        self.page.nodes.append(Node(id, kind, x, y, w, h, label, list(body), parent, **kw))
        return id

    def table(self, id, x, y, w, h, cols, header, rows, parent=None, size=12):
        assert sum(cols) == w, f"table {id}: columns sum {sum(cols)} != width {w}"
        n = Node(id, "panel", x, y, w, h or 10, "", [], parent,
                 table={"cols": cols, "header": header, "rows": rows}, size=size)
        if not h:                       # h=0 means "as tall as the rows need"
            n.h = table_layout(n)[1] + 2
        self.page.nodes.append(n)
        return id

    def edge(self, src, dst, kind="payload", label="", s="r", d="l", via=(), s_at=0.5, d_at=0.5):
        self._n += 1
        self.page.edges.append(Edge(f"e{self._n}", src, dst, kind, label, s, d, [list(p) for p in via], s_at, d_at))

    def chain(self, ids, kind="payload", s="r", d="l", labels=None):
        for i in range(len(ids) - 1):
            self.edge(ids[i], ids[i + 1], kind, (labels or {}).get(i, ""), s, d)


# ------------------------------------------------------------------ geometry helpers
def node_map(page):
    return {n.id: n for n in page.nodes}


def anchor(n: Node, side: str, at: float = 0.5):
    return {"l": (n.x, n.y + n.h * at), "r": (n.x + n.w, n.y + n.h * at),
            "t": (n.x + n.w * at, n.y), "b": (n.x + n.w * at, n.y + n.h)}[side]


def route(a: Node, e: Edge, b: Node):
    """Orthogonal polyline (list of (x, y)) for an edge; `via` points are explicit corners."""
    p0, p1 = anchor(a, e.s, e.s_at), anchor(b, e.d, e.d_at)
    if e.via:
        return [p0] + [tuple(p) for p in e.via] + [p1]
    horiz = lambda s: s in "lr"
    if horiz(e.s) and horiz(e.d):
        if abs(p0[1] - p1[1]) < 0.5:
            return [p0, p1]
        mx = (p0[0] + p1[0]) / 2
        return [p0, (mx, p0[1]), (mx, p1[1]), p1]
    if not horiz(e.s) and not horiz(e.d):
        if abs(p0[0] - p1[0]) < 0.5:
            return [p0, p1]
        my = (p0[1] + p1[1]) / 2
        return [p0, (p0[0], my), (p1[0], my), p1]
    if horiz(e.s):
        return [p0, (p1[0], p0[1]), p1]
    return [p0, (p0[0], p1[1]), p1]


# ------------------------------------------------------------------ text metrics
def _cpl(width, fs, mono):
    return max(4, int((width - 20) / (fs * (0.60 if mono else 0.56))))


def wrap(text, n):
    out, line = [], ""
    for word in text.split(" "):
        while len(word) > n:
            if line:
                out.append(line)
                line = ""
            out.append(word[:n])
            word = word[n:]
        if not line:
            line = word
        elif len(line) + 1 + len(word) <= n:
            line += " " + word
        else:
            out.append(line)
            line = word
    out.append(line)
    return out


def body_lines(n: Node):
    k = _cpl(n.w, n.size, n.mono)
    out = []
    for ln in n.body:
        out.extend(wrap(ln, k) if ln else [""])
    return out


def title_lines(n: Node):
    return wrap(n.label, _cpl(n.w, 16, False)) if n.label else []


def needed_height(n: Node):
    if n.table:
        return table_layout(n)[1]
    if n.kind in CONTAINER_KINDS or n.kind in ("boundary",):
        return 12 + 22 * max(1, len(title_lines(n))) + 8
    if n.kind == "note":
        return int(n.size * 1.35) * len(body_lines(n)) + 2
    return 10 + 21 * len(title_lines(n)) + int(n.size * 1.35) * len(body_lines(n)) + 10


def table_layout(n: Node):
    fs = n.size
    lh = int(fs * 1.35)
    hdr_h = lh + 12
    row_h = []
    for row in n.table["rows"]:
        most = 1
        for ci, cell in enumerate(row):
            k = _cpl(n.table["cols"][ci], fs, False) - 1
            lines = [x for c in ([cell] if isinstance(cell, str) else cell) for x in wrap(c, k)]
            most = max(most, len(lines))
        row_h.append(most * lh + 10)
    return (hdr_h, hdr_h + sum(row_h), row_h)


# ------------------------------------------------------------------ contrast
def _lum(hexc):
    r, g, b = (int(hexc[i:i + 2], 16) / 255 for i in (1, 3, 5))
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def contrast(a, b):
    la, lb = sorted((_lum(a), _lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


# ------------------------------------------------------------------ checks
def check(pages):
    problems, info = [], []
    for pg in pages:
        nm = node_map(pg)
        ids = [n.id for n in pg.nodes]
        if len(ids) != len(set(ids)):
            problems.append(f"[{pg.id}] duplicate node ids")
        for n in pg.nodes:
            if n.x < 0 or n.y < 0 or n.x + n.w > PAGE_W or n.y + n.h > PAGE_H:
                problems.append(f"[{pg.id}] {n.id} outside the {PAGE_W}x{PAGE_H} page: ({n.x},{n.y},{n.w},{n.h})")
            if n.parent:
                p = nm.get(n.parent)
                if p is None:
                    problems.append(f"[{pg.id}] {n.id}: unknown parent {n.parent}")
                elif not (p.x <= n.x and p.y <= n.y and n.x + n.w <= p.x + p.w and n.y + n.h <= p.y + p.h):
                    problems.append(f"[{pg.id}] {n.id} not inside parent {p.id}")
            if n.kind not in KINDS:
                problems.append(f"[{pg.id}] {n.id}: unknown kind {n.kind}")
            need = needed_height(n)
            if n.kind != "note" and need > n.h:
                problems.append(f"[{pg.id}] {n.id}: text needs ~{need}px, box is {n.h}px ({n.label!r})")
            if n.kind == "note" and need > n.h:
                problems.append(f"[{pg.id}] {n.id}: note needs ~{need}px, has {n.h}px")
        # overlap: two nodes overlap only if one contains the other by parentage
        leaf = [n for n in pg.nodes if n.kind not in CONTAINER_KINDS and n.kind != "boundary"]
        for i, a in enumerate(leaf):
            for b in leaf[i + 1:]:
                if a.x < b.x + b.w and b.x < a.x + a.w and a.y < b.y + b.h and b.y < a.y + a.h:
                    problems.append(f"[{pg.id}] overlap: {a.id} and {b.id}")
        for e in pg.edges:
            for end in (e.src, e.dst):
                if end not in nm:
                    problems.append(f"[{pg.id}] edge {e.id} references unknown node {end}")
            if e.kind not in EDGES:
                problems.append(f"[{pg.id}] edge {e.id}: unknown kind {e.kind}")
        info.append(f"{pg.id}: {len(pg.nodes)} nodes, {len(pg.edges)} edges")
    worst = min(contrast(v[2], v[0]) for v in KINDS.values() if v[4] != "text")
    info.append(f"palette: lowest text/fill contrast {worst:.1f}:1 (WCAG AA needs 4.5, AAA 7)")
    return problems, info


# ------------------------------------------------------------------ JSON
def to_json(pages):
    def n2d(n):
        d = asdict(n)
        if not d["table"]:
            del d["table"]
        return {k: v for k, v in d.items() if v not in ("", [], None, False) or k in ("x", "y", "w", "h")}
    return json.dumps({
        "canvas": {"width": PAGE_W, "height": PAGE_H, "grid": 8, "margin": 64},
        "fonts": {"sans": FONT_SANS, "mono": FONT_MONO},
        "node_kinds": {k: {"fill": v[0], "stroke": v[1], "text": v[2], "dashed": v[3], "shape": v[4], "meaning": v[5]}
                       for k, v in KINDS.items()},
        "edge_kinds": {k: {"color": v[0], "width": v[1], "dash": v[2], "meaning": v[3]} for k, v in EDGES.items()},
        "pages": [{"id": p.id, "name": p.name, "title": p.title, "subtitle": p.subtitle, "footer": p.footer, "talk": p.talk,
                   "nodes": [n2d(n) for n in p.nodes], "edges": [asdict(e) for e in p.edges]} for p in pages],
    }, indent=1, ensure_ascii=False)


# ------------------------------------------------------------------ draw.io
def _esc(s):
    return html.escape(s, quote=True)


def _label_html(n: Node):
    if n.table:
        t = n.table
        cols = "".join(f'<col width="{w}"/>' for w in t["cols"])
        head = "".join(f'<th align="left" style="border:1px solid #B0BEC5;background:#ECEFF1;padding:3px 6px">{_esc(h)}</th>' for h in t["header"])
        rows = ""
        for r in t["rows"]:
            cells = ""
            for c in r:
                lines = [c] if isinstance(c, str) else c
                cells += f'<td valign="top" style="border:1px solid #CFD8DC;padding:3px 6px">{"<br/>".join(_esc(x) for x in lines if x)}</td>'
            rows += f"<tr>{cells}</tr>"
        return (f'<table cellspacing="0" style="border-collapse:collapse;width:100%;font-size:{n.size}px;'
                f'font-family:{FONT_SANS}"><colgroup>{cols}</colgroup><tr>{head}</tr>{rows}</table>')
    fam = FONT_MONO if n.mono else FONT_SANS
    out = ""
    if n.label:
        out += f'<b style="font-size:16px">{_esc(n.label)}</b>'
    if n.body:
        out += ("<br/>" if n.label else "") + f'<span style="font-size:{n.size}px;font-family:{fam}">' + "<br/>".join(_esc(b) for b in n.body) + "</span>"
    return out


def _node_style(n: Node):
    fill, stroke, text, dashed, shape, _ = KINDS[n.kind]
    base = f"html=1;whiteSpace=wrap;fontFamily=Inter;fontColor={text};fillColor={fill};strokeColor={stroke};strokeWidth=1.5;"
    if n.table:
        return "html=1;whiteSpace=wrap;overflow=fill;strokeColor=none;fillColor=none;align=left;verticalAlign=top;spacing=0;"
    if shape == "text":
        return "text;html=1;whiteSpace=wrap;overflow=hidden;strokeColor=none;fillColor=none;align=left;verticalAlign=top;spacing=2;" + f"fontColor={text};"
    if shape == "container":
        return base + f"rounded=1;arcSize=2;container=1;collapsible=0;recursiveResize=0;align=left;verticalAlign=top;spacingLeft=14;spacingTop=6;{'dashed=1;strokeWidth=2.5;' if dashed else ''}"
    if shape == "band":
        return base + "rounded=1;arcSize=12;align=center;verticalAlign=middle;strokeWidth=2.5;"
    extra = {"document": "shape=document;boundedLbl=1;", "process": "shape=process;", "rounded": "rounded=1;arcSize=8;"}[shape]
    return base + extra + f"align=left;verticalAlign=top;spacingLeft=10;spacingTop=6;spacingRight=8;{'dashed=1;dashPattern=8 4;' if dashed else ''}"


def to_drawio(pages):
    out = ['<?xml version="1.0" encoding="UTF-8"?>', '<mxfile host="app.diagrams.net" version="24.0.0">']
    for pg in pages:
        nm = node_map(pg)

        def depth(n):
            d = 0
            while n.parent:
                n = nm[n.parent]
                d += 1
            return d

        out.append(f'<diagram id="{_esc(pg.id)}" name="{_esc(pg.name)}">')
        out.append(f'<mxGraphModel dx="1600" dy="900" grid="1" gridSize="8" guides="1" tooltips="1" connect="1" arrows="1" fold="1" '
                   f'page="1" pageScale="1" pageWidth="{PAGE_W}" pageHeight="{PAGE_H}" math="0" shadow="0"><root>')
        out.append('<mxCell id="0"/><mxCell id="1" parent="0"/>')
        # slide title, subtitle, footer are plain text cells
        for cid, text, style, y, h in (
            ("title", f"<b>{_esc(pg.title)}</b>", f"text;html=1;fontSize=32;fontFamily=Inter;fontColor=#0F2A43;align=left;verticalAlign=top;", 28, 44),
            ("subtitle", _esc(pg.subtitle), "text;html=1;fontSize=17;fontFamily=Inter;fontColor=#546E7A;align=left;verticalAlign=top;", 76, 28),
            ("footer", _esc(pg.footer), "text;html=1;fontSize=11;fontFamily=Inter;fontColor=#78909C;align=left;verticalAlign=top;", PAGE_H - 30, 22),
        ):
            out.append(f'<mxCell id="{pg.id}-{cid}" value="{_esc(text)}" style="{style}" vertex="1" parent="1">'
                       f'<mxGeometry x="64" y="{y}" width="{PAGE_W - 128}" height="{h}" as="geometry"/></mxCell>')
        for n in sorted(pg.nodes, key=depth):
            if n.parent:
                p = nm[n.parent]
                gx, gy, par = n.x - p.x, n.y - p.y, f"{pg.id}-{n.parent}"
            else:
                gx, gy, par = n.x, n.y, "1"
            val = _label_html(n)
            out.append(f'<mxCell id="{pg.id}-{n.id}" value="{_esc(val)}" style="{_node_style(n)}" vertex="1" parent="{par}">'
                       f'<mxGeometry x="{gx}" y="{gy}" width="{n.w}" height="{n.h}" as="geometry"/></mxCell>')
        ex = lambda side, at: {"l": (0, at), "r": (1, at), "t": (at, 0), "b": (at, 1)}[side]
        for e in pg.edges:
            color, width, dash, _ = EDGES[e.kind]
            st = (f"edgeStyle=orthogonalEdgeStyle;rounded=1;html=1;endArrow=block;endFill=1;strokeColor={color};strokeWidth={width};"
                  f"exitX={ex(e.s, e.s_at)[0]:.4f};exitY={ex(e.s, e.s_at)[1]:.4f};entryX={ex(e.d, e.d_at)[0]:.4f};entryY={ex(e.d, e.d_at)[1]:.4f};fontSize=12;fontFamily=Inter;"
                  f"labelBackgroundColor=#FFFFFF;fontColor={color};")
            if dash:
                st += "dashed=1;dashPattern=" + dash.replace(",", " ") + ";"
            pts = "".join(f'<mxPoint x="{x}" y="{y}"/>' for x, y in e.via)
            geo = f'<Array as="points">{pts}</Array>' if pts else ""
            out.append(f'<mxCell id="{pg.id}-{e.id}" value="{_esc(e.label)}" style="{st}" edge="1" parent="1" '
                       f'source="{pg.id}-{e.src}" target="{pg.id}-{e.dst}"><mxGeometry relative="1" as="geometry">{geo}</mxGeometry></mxCell>')
        out.append("</root></mxGraphModel></diagram>")
    out.append("</mxfile>")
    return "\n".join(out)


# ------------------------------------------------------------------ SVG
def to_svg(pg: Page):
    nm = node_map(pg)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {PAGE_W} {PAGE_H}" width="{PAGE_W}" height="{PAGE_H}" '
         f'font-family="{_esc(FONT_SANS)}">', "<defs>"]
    for k, (color, *_r) in EDGES.items():
        o.append(f'<marker id="a-{k}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
                 f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>')
    o.append("</defs>")
    o.append(f'<rect width="{PAGE_W}" height="{PAGE_H}" fill="#FFFFFF"/>')
    o.append(f'<text x="64" y="66" font-size="32" font-weight="700" fill="#0F2A43">{_esc(pg.title)}</text>')
    o.append(f'<text x="64" y="98" font-size="17" fill="#546E7A">{_esc(pg.subtitle)}</text>')
    o.append(f'<text x="64" y="{PAGE_H - 14}" font-size="11" fill="#78909C">{_esc(pg.footer)}</text>')

    def depth(n):
        d = 0
        while n.parent:
            n = nm[n.parent]
            d += 1
        return d

    for n in sorted(pg.nodes, key=depth):
        fill, stroke, text, dashed, shape, _ = KINDS[n.kind]
        if n.table:
            hdr_h, _tot, row_h = table_layout(n)
            fs, lh = n.size, int(n.size * 1.35)
            o.append(f'<rect x="{n.x}" y="{n.y}" width="{n.w}" height="{n.h}" fill="#FFFFFF" stroke="#B0BEC5"/>')
            o.append(f'<rect x="{n.x}" y="{n.y}" width="{n.w}" height="{hdr_h}" fill="#ECEFF1" stroke="#B0BEC5"/>')
            cx = n.x
            for ci, w in enumerate(n.table["cols"]):
                o.append(f'<text x="{cx + 6}" y="{n.y + hdr_h - 8}" font-size="{fs}" font-weight="700" fill="#263238">{_esc(n.table["header"][ci])}</text>')
                if ci:
                    o.append(f'<line x1="{cx}" y1="{n.y}" x2="{cx}" y2="{n.y + n.h}" stroke="#CFD8DC"/>')
                cx += w
            ry = n.y + hdr_h
            for ri, row in enumerate(n.table["rows"]):
                o.append(f'<line x1="{n.x}" y1="{ry}" x2="{n.x + n.w}" y2="{ry}" stroke="#CFD8DC"/>')
                cx = n.x
                for ci, cell in enumerate(row):
                    k = _cpl(n.table["cols"][ci], fs, False) - 1
                    lines = [x for c in ([cell] if isinstance(cell, str) else cell) for x in wrap(c, k)]
                    for li, ln in enumerate(l for l in lines if l is not None):
                        o.append(f'<text x="{cx + 6}" y="{ry + 4 + lh * (li + 1) - 4}" font-size="{fs}" fill="#263238">{_esc(ln)}</text>')
                    cx += n.table["cols"][ci]
                ry += row_h[ri]
            continue
        if shape == "text":
            for i, ln in enumerate(body_lines(n)):
                o.append(f'<text x="{n.x + 2}" y="{n.y + 14 + i * int(n.size * 1.35)}" font-size="{n.size}" fill="{text}">{_esc(ln)}</text>')
            continue
        dash = ' stroke-dasharray="8 4"' if dashed else ""
        sw = 2.5 if shape == "band" or n.kind == "zone_gate" else 1.5
        rx = 4 if shape == "container" else 8 if shape != "process" else 3
        o.append(f'<rect x="{n.x}" y="{n.y}" width="{n.w}" height="{n.h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{dash}/>')
        if shape == "process":
            o.append(f'<line x1="{n.x + 10}" y1="{n.y}" x2="{n.x + 10}" y2="{n.y + n.h}" stroke="{stroke}"/>')
            o.append(f'<line x1="{n.x + n.w - 10}" y1="{n.y}" x2="{n.x + n.w - 10}" y2="{n.y + n.h}" stroke="{stroke}"/>')
        tx = n.x + (n.w / 2 if shape == "band" else 14 if shape == "container" else 12 + (10 if shape == "process" else 0))
        anchor_ = "middle" if shape == "band" else "start"
        y = n.y + 22
        for ln in title_lines(n):
            o.append(f'<text x="{tx}" y="{y}" font-size="16" font-weight="700" fill="{text}" text-anchor="{anchor_}">{_esc(ln)}</text>')
            y += 21
        fam = f' font-family="{_esc(FONT_MONO)}"' if n.mono else ""
        for ln in body_lines(n):
            o.append(f'<text x="{tx}" y="{y}" font-size="{n.size}" fill="{text}" text-anchor="{anchor_}"{fam}>{_esc(ln)}</text>')
            y += int(n.size * 1.35)
    for e in pg.edges:
        color, width, dash, _ = EDGES[e.kind]
        pts = route(nm[e.src], e, nm[e.dst])
        d = " ".join(("M" if i == 0 else "L") + f"{x:.0f},{y:.0f}" for i, (x, y) in enumerate(pts))
        da = f' stroke-dasharray="{dash.replace(",", " ")}"' if dash else ""
        o.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"{da} marker-end="url(#a-{e.kind})"/>')
        if e.label:
            seg = max(zip(pts, pts[1:]), key=lambda s: abs(s[0][0] - s[1][0]) + abs(s[0][1] - s[1][1]))
            mx, my = (seg[0][0] + seg[1][0]) / 2, (seg[0][1] + seg[1][1]) / 2
            o.append(f'<text x="{mx:.0f}" y="{my - 5:.0f}" font-size="12" fill="{color}" text-anchor="middle" stroke="#FFFFFF" '
                     f'stroke-width="4" paint-order="stroke">{_esc(e.label)}</text>')
    o.append("</svg>")
    return "\n".join(o)


# ------------------------------------------------------------------ Markdown
def to_markdown(pages):
    md = []
    for pg in pages:
        nm = node_map(pg)
        md.append(f"\n## Page {pg.id[1:]}: {pg.name}\n")
        md.append(f"**Slide title:** {pg.title}  \n**Subtitle:** {pg.subtitle}  \n**Footer:** {pg.footer}\n")
        if pg.talk:
            md.append("**Talk track**\n")
            md.extend(f"- {t}" for t in pg.talk)
            md.append("")
        containers = [n for n in pg.nodes if n.kind in CONTAINER_KINDS or n.kind == "boundary"]
        md.append("**Containers (grouping hierarchy)**\n")
        md.append("| id | kind | x | y | w | h | title | children |")
        md.append("|:--|:--|--:|--:|--:|--:|:--|--:|")
        for c in containers:
            kids = sum(1 for n in pg.nodes if n.parent == c.id)
            md.append(f"| `{c.id}` | {c.kind} | {c.x} | {c.y} | {c.w} | {c.h} | {c.label} | {kids} |")
        md.append("\n**Nodes**\n")
        md.append("| id | kind | parent | x | y | w | h | text (title / body lines) |")
        md.append("|:--|:--|:--|--:|--:|--:|--:|:--|")
        for n in pg.nodes:
            if n in containers:
                continue
            if n.table:
                txt = "TABLE: " + " | ".join(n.table["header"]) + f" ({len(n.table['rows'])} rows)"
            else:
                txt = " / ".join(x for x in ([f"**{n.label}**"] if n.label else []) + [b for b in n.body if b])
            esc_txt = txt.replace("|", "\\|")
            md.append(f"| `{n.id}` | {n.kind} | {n.parent or ''} | {n.x} | {n.y} | {n.w} | {n.h} | {esc_txt} |")
        tables = [n for n in pg.nodes if n.table]
        for t in tables:
            md.append(f"\n**Table `{t.id}`** (columns {t.table['cols']} px)\n")
            md.append("| " + " | ".join(t.table["header"]) + " |")
            md.append("|" + "|".join(":--" for _ in t.table["header"]) + "|")
            for r in t.table["rows"]:
                md.append("| " + " | ".join(("<br>".join(x for x in ([c] if isinstance(c, str) else c) if x)).replace("|", "\\|") for c in r) + " |")
        md.append("\n**Edges (directed)**\n")
        md.append("| id | from | exit | to | entry | kind | label | via |")
        md.append("|:--|:--|:-:|:--|:-:|:--|:--|:--|")
        for e in pg.edges:
            via = " > ".join(f"({x},{y})" for x, y in e.via)
            md.append(f"| {e.id} | `{e.src}` | {e.s} | `{e.dst}` | {e.d} | {e.kind} | {e.label} | {via} |")
    return "\n".join(md) + "\n"
