"""poster/poster.md -> poster/out/poster.pdf (+ PNG preview). Run: python poster/build.py

Same idea as a print pipeline: the manuscript is canonical and read-only, every coordinate comes
from SPEC (mm, page origin top-left), and the build ends with machine checks. Looking at the PNG
is still required; the checks only catch what a machine can see.

Checks
1. every paragraph, bullet and figure carries an {E..} tag, and each tagged Entry exists in
   notes/lab_notebook.md (nothing goes on the poster without its entry);
2. every figure file exists;
3. no column runs past the bottom margin (no silent reflow: move "@column" instead);
4. no text below SPEC["min_pt"].
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# Word-by-word placement needs text widths that scale linearly with dpi; hinting breaks that and
# makes long words overrun the next one in a preview rendered at another dpi.
matplotlib.rcParams.update({"text.hinting": "none", "text.hinting_factor": 1})
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "poster" / "poster.md"
OUT = ROOT / "poster" / "out"
NOTEBOOK = ROOT / "notes" / "lab_notebook.md"

# ponytail: A0 portrait is a default, not a KPS rule. Check the conference poster size and edit here.
SPEC = dict(
    W=841, H=1189, margin=24, gap=18, cols=3,
    header_h=118,
    title_pt=66, authors_pt=34, head_pt=38, body_pt=27, caption_pt=21, min_pt=20,
    leading=1.32, para_gap=7, panel_gap=16, fig_gap=5,
    ink="#1b1b1b", ink2="#55524c", accent="#2a5d9f", band="#eef2f8", rule="#c9d3e3",
    font="DejaVu Sans",
)
PT_MM = 25.4 / 72


# ───────────── manuscript ─────────────

TAG = re.compile(r"\{E(\d+)\}")


def parse(text):
    meta, body = {}, text
    if text.startswith("---"):
        _, fm, body = text.split("---", 2)
        for line in fm.strip().splitlines():
            k, v = line.split(":", 1)
            meta[k.strip()] = v.strip()
    body = re.sub(r"<!--.*?-->", "", body, flags=re.S)
    panels, cur = [], None
    for block in re.split(r"\n\s*\n", body.strip()):
        block = block.strip()
        if block == "@column":
            panels.append("@column")
            continue
        if block.startswith("## "):
            num, title = [s.strip() for s in block[3:].split("|", 1)]
            cur = dict(num=num, title=title, items=[])
            panels.append(cur)
            continue
        m = re.match(r"!\[(.*?)\]\((.*?)\)(.*)", block, re.S)
        if m:
            cur["items"].append(dict(kind="fig", caption=m[1], path=m[2], tags=TAG.findall(m[3])))
            continue
        if block.startswith("- "):
            for b in re.split(r"\n- ", block[2:]):
                b = " ".join(b.split())
                cur["items"].append(dict(kind="bullet", text=TAG.sub("", b).strip(), tags=TAG.findall(b)))
            continue
        t = " ".join(block.split())
        cur["items"].append(dict(kind="para", text=TAG.sub("", t).strip(), tags=TAG.findall(t)))
    return meta, panels


# ───────────── typesetting ─────────────

TOKEN = re.compile(r"\$[^$]+\$|\*\*[^*]+\*\*|[^\s$*]+")


class Page:
    def __init__(self):
        s = SPEC
        self.fig = plt.figure(figsize=(s["W"] / 25.4, s["H"] / 25.4), dpi=100)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, s["W"])
        self.ax.set_ylim(s["H"], 0)
        self.ax.axis("off")
        self.r = self.fig.canvas.get_renderer()
        self._w = {}

    def width(self, tok, pt, bold=False):
        key = (tok, pt, bold)
        if key not in self._w:
            t = self.ax.text(0, 0, tok, fontsize=pt, fontweight="bold" if bold else "normal", family=SPEC["font"])
            self._w[key] = t.get_window_extent(self.r).width / self.fig.dpi * 25.4
            t.remove()
        return self._w[key]

    def flow(self, text, x, y, w, pt, color=None, bold=False, indent=0.0):
        """Word-by-word placement with wrapping; $..$ and **..** stay atomic. Returns the new y (mm)."""
        color = color or SPEC["ink"]
        lh = pt * PT_MM * SPEC["leading"]
        sp = self.width("a a", pt) - 2 * self.width("a", pt)
        line, lw, y = [], 0.0, y + pt * PT_MM          # baseline of the first line
        prev = 0
        toks = []
        for m in TOKEN.finditer(text):
            glue = m.start() == prev and toks
            prev = m.end()
            s, b = m[0], bold
            if s.startswith("**"):
                s, b = s[2:-2], True
            toks.append((s, b, bool(glue)))

        def emit(line, y):
            cx = x + indent
            for s, b, g, tw in line:
                if not g and cx > x + indent:
                    cx += sp
                self.ax.text(cx, y, s, fontsize=pt, color=color, family=SPEC["font"],
                             fontweight="bold" if b else "normal", va="baseline")
                cx += tw
        for s, b, g in toks:
            tw = self.width(s, pt, b)
            add = tw + (0 if g or not line else sp)
            if line and lw + add > w - indent:
                emit(line, y)
                y += lh
                line, lw, add = [], 0.0, tw
                g = False
            line.append((s, b, g, tw))
            lw += add
        if line:
            emit(line, y)
        return y + lh - pt * PT_MM


def build():
    s = SPEC
    meta, panels = parse(MS.read_text(encoding="utf8"))
    known = set(re.findall(r"^## .*?Entry (\d+)\b", NOTEBOOK.read_text(encoding="utf8"), flags=re.M))
    problems, used = [], set()
    P = Page()
    ax = P.ax

    # header band
    ax.add_patch(plt.Rectangle((0, 0), s["W"], s["header_h"], color=s["band"], lw=0))
    ax.add_patch(plt.Rectangle((0, s["header_h"] - 2.5), s["W"], 2.5, color=s["accent"], lw=0))
    y = P.flow(meta.get("title", ""), s["margin"], s["margin"] - 4, s["W"] - 2 * s["margin"], s["title_pt"], bold=True)
    y = P.flow(meta.get("authors", ""), s["margin"], y + 5, s["W"] - 2 * s["margin"], s["authors_pt"], s["ink2"])
    P.flow(" · ".join(v for k, v in meta.items() if k in ("affiliation", "contact")), s["margin"], y + 2,
           s["W"] - 2 * s["margin"], s["authors_pt"] - 6, s["ink2"])

    cw = (s["W"] - 2 * s["margin"] - (s["cols"] - 1) * s["gap"]) / s["cols"]
    col, top = 0, s["header_h"] + s["panel_gap"]
    y = top
    bottom = s["H"] - s["margin"]
    fill = []

    for p in panels:
        if p == "@column":
            fill.append(y)
            col, y = col + 1, top
            continue
        x = s["margin"] + col * (cw + s["gap"])
        ax.text(x, y + s["head_pt"] * PT_MM, p["num"], fontsize=s["head_pt"], color=s["accent"],
                fontweight="bold", family=s["font"], va="baseline")
        nx = x + P.width(p["num"], s["head_pt"], True) + 4
        y = P.flow(p["title"], nx, y, cw - (nx - x), s["head_pt"], bold=True)
        ax.plot([x, x + cw], [y + 2.5, y + 2.5], color=s["rule"], lw=2)
        y += 2.5 + s["para_gap"]
        for it in p["items"]:
            if not it["tags"]:
                problems.append(f"panel {p['num']}: untagged {it['kind']}: {it.get('text', it.get('caption'))[:50]!r}")
            for t in it["tags"]:
                used.add(int(t))
                if t not in known:
                    problems.append(f"panel {p['num']}: tag E{t} is not an Entry in the notebook")
            if it["kind"] == "fig":
                f = ROOT / it["path"]
                if not f.exists():
                    problems.append(f"panel {p['num']}: missing figure {it['path']}")
                    continue
                img = mpimg.imread(f)
                h = cw * img.shape[0] / img.shape[1]
                fa = P.fig.add_axes([x / s["W"], 1 - (y + h) / s["H"], cw / s["W"], h / s["H"]])
                fa.imshow(img, interpolation="lanczos")
                fa.axis("off")
                y += h + s["fig_gap"]
                y = P.flow(it["caption"], x, y, cw, s["caption_pt"], s["ink2"]) + s["para_gap"]
            elif it["kind"] == "bullet":
                ax.text(x, y + s["body_pt"] * PT_MM, "•", fontsize=s["body_pt"], color=s["accent"],
                        family=s["font"], va="baseline")
                y = P.flow(it["text"], x, y, cw, s["body_pt"], indent=7) + s["para_gap"]
            else:
                y = P.flow(it["text"], x, y, cw, s["body_pt"]) + s["para_gap"]
        y += s["panel_gap"]
        if y - s["panel_gap"] > bottom:
            problems.append(f"column {col + 1} overflows the bottom margin by {y - s['panel_gap'] - bottom:.0f} mm "
                            f"(at panel {p['num']})")
    fill.append(y)
    for k in ("body_pt", "caption_pt", "head_pt"):
        if s[k] < s["min_pt"]:
            problems.append(f"{k} = {s[k]} pt < min {s['min_pt']} pt")

    OUT.mkdir(exist_ok=True)
    P.fig.savefig(OUT / "poster.pdf")
    P.fig.savefig(OUT / "poster.png", dpi=60)
    avail = bottom - top
    print("columns used: " + "  ".join(f"{i + 1}: {100 * (f - s['panel_gap'] - top) / avail:.0f}%"
                                       for i, f in enumerate(fill)))
    print("entries cited: " + ", ".join(f"E{e}" for e in sorted(used)))
    print(f"wrote {OUT / 'poster.pdf'} and poster.png")
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for m in problems:
            print("  ! " + m)
        return 1
    print("checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(build())
