"""Prose for graphic-spec.md; the palette and edge tables are generated from graphics.py."""

HEADER = """# pto-core: editable graphic specification

Five 16:9 slides on a 1920 x 1080 canvas. Everything here is generated from one model,
`src/pages.py`, so the machine files and this document cannot disagree:

| File | Use it to |
|:--|:--|
| `pto-core-diagrams.drawio` | open all five pages in draw.io / diagrams.net (File > Open from > Device). Containers are real parent shapes, so moving one moves its children; edges stay connected to their nodes. |
| `svg/*.svg` | upload to Canva (Uploads, drag the five files), or place in PowerPoint, Keynote or Figma. Text stays live text. |
| `diagram-spec.json` | rebuild in Excalidraw or any scripting tool: `pages[].nodes` are rectangles (x, y, w, h, kind, text) and `pages[].edges` are arrows (from, to, exit side, entry side, via points). |

Regenerate after editing `src/pages.py`:

```bash
python3 docs/diagrams/src/build_graphics.py docs/diagrams
```

The build fails, and writes nothing, if a node leaves the page, overlaps a sibling or sits outside its
parent, if an edge names a missing node, or if a label needs more height than its box has.

## Canvas, grid and type

| Item | Value |
|:--|:--|
| Slide | 1920 x 1080 px, 16:9; outer margin 64 px; snap to an 8 px grid |
| Title band | y 28-104: title 32 px bold `#0F2A43`, subtitle 17 px `#546E7A` (one sentence) |
| Footer | y 1050, 11 px `#78909C`: source files and the commit the claims were read at |
| Text face | Inter (Helvetica or Arial as fallback); code identifiers may use JetBrains Mono |
| Sizes | container and node titles 16 px bold; body 13 px; table cells 12 px; edge labels 12 px |
| Node anatomy | title top-left, 10 px inset; body lines below, no bullets; 1.5 px stroke; radius 8 px |
| Container anatomy | title top-left in the container's own colour; children inset at least 24 px |
"""

SHAPES = """
## Shape language

| Shape | Meaning |
|:--|:--|
| rounded rectangle | a stage, component or gate |
| folded-corner document | an artifact, file format or stream |
| rectangle with two vertical bars | a concurrency primitive: ring, pool, lock |
| wide horizontal band | the process boundary |
| dashed outline | outside this repository: its contract is defined here, its code is not |
| large container | grouping only; carries no behaviour |

## Accessibility

Meaning never rides on colour alone: each kind also differs by outline (solid, dashed) or shape (document,
double bar, band), and every engine node states "C++20" or "CLI" in its own text. Text-on-fill contrast is
checked at build time; the lowest pair in this palette is __WORST__:1 (WCAG AAA needs 7:1). The smallest
text is 11 px and is used only for footers.
"""
