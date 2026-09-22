"""A tiny ASCII/Unicode box-drawing canvas.

Boxes and connectors are stored as per-cell direction masks and rendered to the
right junction glyph, so corners, tees and crossings are correct by construction
and a box border can never be misaligned. Text is overlaid last.

Only box-drawing (U+2500..257F) and ASCII are emitted, plus the four arrowheads
the brief names (U+25BA, U+25C4, U+25B2, U+25BC).
"""
from __future__ import annotations

N, E, S, W = "N", "E", "S", "W"
DIRS = (N, E, S, W)

# (N, E, S, W) each 0 | 'l' (light) | 'd' (double)  ->  glyph
_G: dict[tuple, str] = {}


def _reg(glyph: str, n=0, e=0, s=0, w=0):
    _G[(n, e, s, w)] = glyph


# light
for g, k in {
    "│": ("l", 0, "l", 0), "─": (0, "l", 0, "l"),
    "┌": (0, "l", "l", 0), "┐": (0, 0, "l", "l"), "└": ("l", "l", 0, 0), "┘": ("l", 0, 0, "l"),
    "├": ("l", "l", "l", 0), "┤": ("l", 0, "l", "l"), "┬": (0, "l", "l", "l"), "┴": ("l", "l", 0, "l"),
    "┼": ("l", "l", "l", "l"),
}.items():
    _reg(g, *k)
# double
for g, k in {
    "║": ("d", 0, "d", 0), "═": (0, "d", 0, "d"),
    "╔": (0, "d", "d", 0), "╗": (0, 0, "d", "d"), "╚": ("d", "d", 0, 0), "╝": ("d", 0, 0, "d"),
    "╠": ("d", "d", "d", 0), "╣": ("d", 0, "d", "d"), "╦": (0, "d", "d", "d"), "╩": ("d", "d", 0, "d"),
    "╬": ("d", "d", "d", "d"),
}.items():
    _reg(g, *k)
# mixed single/double
for g, k in {
    "╒": (0, "d", "l", 0), "╓": (0, "l", "d", 0), "╕": (0, 0, "l", "d"), "╖": (0, 0, "d", "l"),
    "╘": ("l", "d", 0, 0), "╙": ("d", "l", 0, 0), "╛": ("l", 0, 0, "d"), "╜": ("d", 0, 0, "l"),
    "╞": ("l", "d", "l", 0), "╟": ("d", "l", "d", 0), "╡": ("l", 0, "l", "d"), "╢": ("d", 0, "d", "l"),
    "╤": (0, "d", "l", "d"), "╥": (0, "l", "d", "l"), "╧": ("l", "d", 0, "d"), "╨": ("d", "l", 0, "l"),
    "╪": ("l", "d", "l", "d"), "╫": ("d", "l", "d", "l"),
}.items():
    _reg(g, *k)

_DASH = {"─": "┄", "│": "┆"}
_HEAD = {"E": "►", "W": "◄", "N": "▲", "S": "▼"}
ALLOWED_EXTRA = set("►◄▲▼")


class Canvas:
    problems: list | None = None   # set to [] to collect instead of raising

    def _fail(self, msg: str):
        if self.problems is None:
            raise AssertionError(msg)
        self.problems.append(msg)

    def __init__(self, width: int, height: int):
        self.w, self.h = width, height
        self.mask = [[dict() for _ in range(width)] for _ in range(height)]
        self.dashed = [[False] * width for _ in range(height)]
        self.txt: dict[tuple[int, int], str] = {}

    # ------------------------------------------------------------ primitives
    def _add(self, x: int, y: int, d: str, style: str):
        assert 0 <= x < self.w and 0 <= y < self.h, f"({x},{y}) outside {self.w}x{self.h}"
        cur = self.mask[y][x].get(d)
        # double beats light when two styles claim one arm of one cell
        if cur is None or (cur == "l" and style == "d"):
            self.mask[y][x][d] = style

    def hline(self, x1: int, x2: int, y: int, style: str = "l", dashed: bool = False):
        if x1 > x2:
            x1, x2 = x2, x1
        for x in range(x1, x2 + 1):
            if x > x1:
                self._add(x, y, W, style)
            if x < x2:
                self._add(x, y, E, style)
            if dashed:
                self.dashed[y][x] = True

    def vline(self, x: int, y1: int, y2: int, style: str = "l", dashed: bool = False):
        if y1 > y2:
            y1, y2 = y2, y1
        for y in range(y1, y2 + 1):
            if y > y1:
                self._add(x, y, N, style)
            if y < y2:
                self._add(x, y, S, style)
            if dashed:
                self.dashed[y][x] = True

    def text(self, x: int, y: int, s: str):
        for i, ch in enumerate(s):
            if ch != "\0":
                if not (0 <= x + i < self.w and 0 <= y < self.h):
                    self._fail(f"text at ({x},{y}) len {len(s)} overflows width {self.w}: {s[:60]!r}")
                    return
                self.txt[(x + i, y)] = ch

    def head(self, x: int, y: int, d: str):
        self.txt[(x, y)] = _HEAD[d]

    # --------------------------------------------------------------- shapes
    def box(self, x: int, y: int, w: int, h: int, title: str = "", lines=(), style: str = "l",
            dashed: bool = False, align: str = "left", pad: int = 2):
        x2, y2 = x + w - 1, y + h - 1
        self.hline(x, x2, y, style, dashed)
        self.hline(x, x2, y2, style, dashed)
        self.vline(x, y, y2, style, dashed)
        self.vline(x2, y, y2, style, dashed)
        if title:
            t = f" {title} "
            if len(t) > w - 3:
                self._fail(f"title {title!r} ({len(t)}) too long for width {w}")
            self.text(x + 2, y, t[: max(0, w - 3)])
        inner = w - 2 * pad
        for i, ln in enumerate(lines):
            if len(ln) > inner:
                self._fail(f"box {title[:28]!r}: line ({len(ln)}>{inner}) {ln[:50]!r}")
                ln = ln[:inner]
            if not (y + 1 + i < y2):
                self._fail(f"box {title[:28]!r} has too many lines ({len(lines)}) for height {h}")
                continue
            if align == "center":
                ln = ln.center(inner)
            self.text(x + pad, y + 1 + i, ln)

    def path(self, pts, style: str = "l", head: str | None = "end", dashed: bool = False):
        """Orthogonal polyline through pts. `head` = 'end' | 'start' | 'both' | None."""
        for (xa, ya), (xb, yb) in zip(pts, pts[1:]):
            assert xa == xb or ya == yb, f"non-orthogonal segment {(xa, ya)}->{(xb, yb)}"
            if ya == yb:
                self.hline(xa, xb, ya, style, dashed)
            else:
                self.vline(xa, ya, yb, style, dashed)
        if head in ("end", "both"):
            (xa, ya), (xb, yb) = pts[-2], pts[-1]
            d = E if xb > xa else W if xb < xa else S if yb > ya else N
            self.head(xb, yb, d)
        if head in ("start", "both"):
            (xa, ya), (xb, yb) = pts[1], pts[0]
            d = E if xb > xa else W if xb < xa else S if yb > ya else N
            self.head(xb, yb, d)

    # -------------------------------------------------------------- render
    def _glyph(self, x: int, y: int) -> str:
        m = self.mask[y][x]
        if not m:
            return " "
        key = tuple(m.get(d, 0) for d in DIRS)
        g = _G.get(key)
        if g is None:  # unsupported weight mix: promote to the heavier style
            style = "d" if "d" in m.values() else "l"
            g = _G.get(tuple(style if m.get(d) else 0 for d in DIRS))
            if g is None:  # single arm: draw the straight glyph
                g = "═" if style == "d" and (E in m or W in m) else "║" if style == "d" else "─" if (E in m or W in m) else "│"
        if self.dashed[y][x] and g in _DASH:
            g = _DASH[g]
        return g

    def render(self) -> list[str]:
        rows = []
        for y in range(self.h):
            row = []
            for x in range(self.w):
                row.append(self.txt.get((x, y)) or self._glyph(x, y))
            rows.append("".join(row).rstrip())
        while rows and not rows[-1]:
            rows.pop()
        return rows


def check_render(rows: list[str], max_width: int = 100) -> list[str]:
    """Lint: width and glyph repertoire."""
    problems = []
    for i, r in enumerate(rows, 1):
        if len(r) > max_width:
            problems.append(f"line {i}: width {len(r)} > {max_width}")
        for ch in r:
            o = ord(ch)
            if ch == " " or (32 < o < 127) or 0x2500 <= o <= 0x257F or ch in ALLOWED_EXTRA:
                continue
            problems.append(f"line {i}: disallowed glyph {ch!r} U+{o:04X}")
            break
    return problems


if __name__ == "__main__":  # smoke test
    c = Canvas(60, 14)
    c.box(0, 0, 22, 5, "reader", ["memcpy 4 members", "push lane k mod D"])
    c.box(30, 0, 22, 5, "inflater", ["libdeflate", "1 thread each"], dashed=True)
    c.path([(22, 2), (29, 2)], style="d")
    c.box(0, 8, 20, 4, "sink", ["fd"])
    c.path([(10, 5), (10, 7)], style="l")
    c.path([(41, 5), (41, 6), (10, 6)], style="l", head=None)  # joins the vertical -> tee
    print("\n".join(c.render()))
    print(check_render(c.render()))
