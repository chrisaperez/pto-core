#!/usr/bin/env python3
"""Build the graphic specification: JSON spec, multi-page draw.io file, per-page SVG, Markdown tables.

    python3 build_graphics.py OUT_DIR [--check-only]

Writes OUT_DIR/diagram-spec.json, OUT_DIR/pto-core-diagrams.drawio, OUT_DIR/svg/<page>.svg and
OUT_DIR/graphic-spec.md. Exits non-zero if any geometry or text-fit check fails.
"""
import sys
import xml.dom.minidom
from pathlib import Path

import graphics as g
from pages import PAGES
from design_text import HEADER, SHAPES


def design_md():
    rows = ["\n## Colour language\n", "| kind | fill | stroke | text | outline | contrast | meaning |", "|:--|:--|:--|:--|:--|--:|:--|"]
    worst = 99.0
    for k, (fill, stroke, text, dashed, shape, meaning) in g.KINDS.items():
        if shape == "text":
            continue
        c = g.contrast(text, fill)
        worst = min(worst, c)
        rows.append(f"| `{k}` | `{fill}` | `{stroke}` | `{text}` | {'dashed' if dashed else 'solid'}, {shape} | {c:.1f}:1 | {meaning} |")
    rows += ["\n## Edge language\n", "| kind | colour | width px | dash | meaning |", "|:--|:--|--:|:--|:--|"]
    for k, (color, width, dash, meaning) in g.EDGES.items():
        rows.append(f"| `{k}` | `{color}` | {width} | {dash or 'solid'} | {meaning} |")
    rows.append("\nArrowheads are filled triangles, 7 px. Edge labels are 12 px in the edge colour on a white halo, used only where the gap is wider than the label.")
    return HEADER + "\n".join(rows) + "\n" + SHAPES.replace("__WORST__", f"{worst:.1f}")


def main(argv):
    out = Path(argv[1])
    pages = [p() for p in PAGES]
    problems, info = g.check(pages)
    print("\n".join(info))
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for p in problems:
            print("  -", p)
    if "--check-only" in argv:
        return 1 if problems else 0
    if problems:
        return 1
    (out / "svg").mkdir(parents=True, exist_ok=True)
    (out / "diagram-spec.json").write_text(g.to_json(pages), encoding="utf-8")
    drawio = g.to_drawio(pages)
    xml.dom.minidom.parseString(drawio.encode("utf-8"))            # must be well-formed XML
    (out / "pto-core-diagrams.drawio").write_text(drawio, encoding="utf-8")
    for pg in pages:
        svg = g.to_svg(pg)
        xml.dom.minidom.parseString(svg.encode("utf-8"))
        (out / "svg" / f"{pg.id}-{pg.name.split(' ', 1)[1].lower().replace(' ', '-')}.svg").write_text(svg, encoding="utf-8")
    (out / "graphic-spec.md").write_text(design_md() + g.to_markdown(pages), encoding="utf-8")
    print("wrote spec, drawio, 5 svg and graphic-spec.md to", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
