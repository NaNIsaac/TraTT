"""poster/poster.md -> poster/out/poster.pdf (+ PNG preview). Run: python poster/build.py

Same idea as a print pipeline: the manuscript is canonical and read-only, every coordinate comes
from SPEC (mm, page origin top-left), and the build ends with machine checks. Looking at the PNG
is still required; the checks only catch what a machine can see.

Look: the group's previous A0 poster (ICAMD 2025) — navy header with logo and a three-colour
stripe, light page, white cards, "■ 01. Title" section heads, two columns, yellow callouts.

Manuscript syntax
  front matter: title, authors, affiliation, contact, footer
  "## Title"      a panel (numbered automatically)        "@wide"   next panel spans both columns
  "@column"       continue in the next column             "> text"  yellow callout
  "- text"        bullet                                  "![caption](path)" figure (width = card)
  {E12}           notebook entry tag (stripped from print, required on every item)

Usage: python poster/build.py [preview_dpi]   (default 60; use ~150 to inspect details)

Checks
1. every paragraph, bullet, callout and figure carries an {E..} tag, and each tagged Entry exists in
   notes/lab_notebook.md (nothing goes on the poster without its entry);
2. every figure file exists;
3. no column runs into the footer (no silent reflow: move "@column" instead);
4. no text below SPEC["min_pt"].

Layout runs twice: the first, unsaved pass measures each column; the second spreads the spare
height evenly between the panels of that column, so both columns end level with the footer.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# Word-by-word placement needs text widths that scale linearly with dpi; hinting breaks that and
# makes long words overrun the next one in a preview rendered at another dpi.
matplotlib.rcParams.update({"text.hinting": "none", "text.hinting_factor": 1, "mathtext.fontset": "custom",
                            "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                            "mathtext.sf": "Arial",
                            "pdf.fonttype": 42})   # TrueType in the PDF: editable, accepted by print shops
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "poster" / "poster.md"
OUT = ROOT / "poster" / "out"
ASSETS = ROOT / "poster" / "assets"
NOTEBOOK = ROOT / "notes" / "lab_notebook.md"

# A0 portrait, same as the ICAMD 2025 poster (841 x 1189 mm). Check the KPS rule before printing.
SPEC = dict(
    W=841, H=1189, margin=16, gap=14, cols=2, pad=12,
    header_h=130, stripe_h=5, footer_h=56,
    title_pt=64, authors_pt=34, head_pt=44, body_pt=28, callout_pt=32, caption_pt=22, min_pt=20,
    leading=1.30, para_gap=7, panel_gap=10, fig_gap=4,
    navy="#0f2a3a", page="#efefef", card="#ffffff", head="#2f5597", ink="#1b1b1b", ink2="#55524c",
    callout="#fff2cc", stripe=("#3a3a3a", "#e04b3a", "#f08a3c"), rule="#d5dbe5",
    font="Arial",
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
    panels, cur, wide = [], None, False
    for block in re.split(r"\n\s*\n", body.strip()):
        block = block.strip()
        if block in ("@column", "@wide"):
            if block == "@column":
                panels.append("@column")
            wide = block == "@wide"
            continue
        if block.startswith("## "):
            cur = dict(title=block[3:].strip(), items=[], wide=wide)
            wide = False
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
        kind = "callout" if block.startswith(">") else "para"
        t = " ".join(re.sub(r"^>\s?", "", block, flags=re.M).split())
        cur["items"].append(dict(kind=kind, text=TAG.sub("", t).strip(), tags=TAG.findall(t)))
    return meta, panels


# ───────────── typesetting ─────────────

TOKEN = re.compile(r"\$[^$]+\$|\*\*[^*]+\*\*|[^\s$*]+|\*")   # a lone * (corresponding author) is a word


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

    def flow(self, text, x, y, w, pt, color=None, bold=False, indent=0.0, align="left", z=3, draw=True):
        """Word-by-word placement with wrapping; $..$ and **..** stay atomic. Returns the new y (mm)."""
        color = color or SPEC["ink"]
        lh = pt * PT_MM * SPEC["leading"]
        sp = self.width("a a", pt) - 2 * self.width("a", pt)
        y = y + pt * PT_MM                                  # baseline of the first line
        toks, prev = [], 0
        for m in TOKEN.finditer(text):
            glue = m.start() == prev and bool(toks)
            prev = m.end()
            s, b = m[0], bold
            if s.startswith("**"):
                s, b = s[2:-2], True
            toks.append((s, b, glue))
        lines, line, lw = [], [], 0.0
        for s, b, g in toks:
            tw = self.width(s, pt, b)
            add = tw + (0 if g or not line else sp)
            if line and lw + add > w - indent:
                k = len(line)                                # glued run ending the line moves along
                while g and k > 1 and line[k - 1][2]:
                    k -= 1
                if g:
                    k -= 1
                carry = line[k:] if k > 0 else []
                kept = line[:k] if k > 0 else line
                lines.append((kept, sum(t[3] for t in kept) + sp * sum(1 for t in kept[1:] if not t[2])))
                line = [(cs, cb, i > 0 and cg, ctw) for i, (cs, cb, cg, ctw) in enumerate(carry)]
                lw = sum(t[3] for t in line) + sp * sum(1 for t in line[1:] if not t[2])
                add = tw if not line else tw + (0 if g else sp)
            line.append((s, b, g if line else False, tw))
            lw += add
        if line:
            lines.append((line, lw))
        for line, lw in lines:
            cx = x + indent + ((w - indent - lw) / 2 if align == "center" else 0)
            first = True
            for s, b, g, tw in line:
                if not g and not first:
                    cx += sp
                if draw:
                    self.ax.text(cx, y, s, fontsize=pt, color=color, family=SPEC["font"], zorder=z,
                                 fontweight="bold" if b else "normal", va="baseline")
                cx += tw
                first = False
            y += lh
        return y - lh + pt * PT_MM * 0.35 if lines else y - pt * PT_MM

    def image(self, path, x, y, w, h=None, z=2):
        """Place an image at (x, y) mm, width w; height from the aspect ratio unless given."""
        img = mpimg.imread(path)
        h = h or w * img.shape[0] / img.shape[1]
        a = self.fig.add_axes([x / SPEC["W"], 1 - (y + h) / SPEC["H"], w / SPEC["W"], h / SPEC["H"]], zorder=z)
        a.imshow(img, interpolation="lanczos")
        a.axis("off")
        return h


def build(preview_dpi=60, extra=None, save=True):
    """One layout pass. extra[c] is added after each panel of column c except the last, so the
    columns end level with the footer (filled in by a first, unsaved pass; see __main__)."""
    s = SPEC
    extra = extra or [0.0] * SPEC["cols"]
    meta, panels = parse(MS.read_text(encoding="utf8"))
    known = set(re.findall(r"^## .*?Entry (\d+)\b", NOTEBOOK.read_text(encoding="utf8"), flags=re.M))
    problems, used = [], set()
    P = Page()
    ax = P.ax
    W, H, M = s["W"], s["H"], s["margin"]

    # page, header, stripe
    ax.add_patch(plt.Rectangle((0, 0), W, H, color=s["page"], lw=0, zorder=0))
    ax.add_patch(plt.Rectangle((0, 0), W, s["header_h"], color=s["navy"], lw=0, zorder=1))
    seg = (0.0, 0.58, 0.60, 0.84, 0.86, 1.0)
    for i, c in enumerate(s["stripe"]):
        ax.add_patch(plt.Rectangle((seg[2 * i] * W, s["header_h"] + 2), (seg[2 * i + 1] - seg[2 * i]) * W,
                                   s["stripe_h"], color=c, lw=0, zorder=1))
    logo_w = 0
    if (ASSETS / "logo_ihs.png").exists():
        logo_w = 80
        P.image(ASSETS / "logo_ihs.png", W - M - logo_w, 22, logo_w * 0.7, h=logo_w * 1.0)
    tw = W - 2 * M - 2 * logo_w
    y = P.flow(meta.get("title", ""), M + logo_w, 20, tw, s["title_pt"], "white", bold=True, align="center")
    y = P.flow(meta.get("authors", ""), M + logo_w, y + 8, tw, s["authors_pt"], "white", align="center")
    P.flow(meta.get("affiliation", ""), M + logo_w, y + 4, tw, s["authors_pt"] - 6, "#d8e2ea", align="center")

    # footer
    fy = H - s["footer_h"]
    ax.plot([M, W - M], [fy, fy], color=s["rule"], lw=2, zorder=1)
    lw_ = 0
    if (ASSETS / "logo_cupt.png").exists():
        lw_ = 150
        P.image(ASSETS / "logo_cupt.png", W - M - lw_, fy + 8, lw_)
    P.flow(" · ".join(v for k, v in meta.items() if k in ("footer", "contact")), M, fy + 10,
           W - 2 * M - lw_ - 10, s["caption_pt"], s["ink2"])

    cw = (W - 2 * M - (s["cols"] - 1) * s["gap"]) / s["cols"]
    top = s["header_h"] + s["stripe_h"] + 2 + s["panel_gap"]
    bottom = fy - s["panel_gap"]
    ys = [top] * s["cols"]
    count = [0] * s["cols"]
    col = 0
    num = 0
    for p in panels:
        if p == "@column":
            col += 1
            continue
        num += 1
        if p["wide"]:
            x, w, y0 = M, W - 2 * M, max(ys)
        else:
            x, w, y0 = M + col * (cw + s["gap"]), cw, ys[col]
        ix, iw = x + s["pad"], w - 2 * s["pad"]
        y = y0 + s["pad"] * 0.6
        # section head: ■ 01. Title
        cap = s["head_pt"] * PT_MM * 0.72                     # Arial cap height
        base = y + s["head_pt"] * PT_MM
        ax.add_patch(plt.Rectangle((x, base - cap), cap, cap, color=s["ink"], lw=0, zorder=3))
        hx = x + cap + 5
        y = P.flow(f"{num:02d}. {p['title']}", hx, y, w - (hx - x), s["head_pt"], s["head"], bold=True)
        y += 3
        card_top = y
        y += s["pad"]
        # A wide panel with a figure is laid out side by side: text left, figure right (half width each).
        side = p["wide"] and any(it["kind"] == "fig" for it in p["items"])
        if side:
            fx, ix, iw = ix + iw / 2 + s["gap"] / 2, ix, iw / 2 - s["gap"] / 2
            fh = sum(iw * (lambda a: a.shape[0] / a.shape[1])(mpimg.imread(ROOT / it["path"]))
                     for it in p["items"] if it["kind"] == "fig" and (ROOT / it["path"]).exists())
            th = text_height(P, p["items"], ix, iw) - s["para_gap"]
            fy = y + max(0.0, (th - fh) / 2)
            y += max(0.0, (fh - th) / 2)
        for it in p["items"]:
            if not it["tags"]:
                problems.append(f"panel {num:02d}: untagged {it['kind']}: {it.get('text', it.get('caption'))[:50]!r}")
            for t in it["tags"]:
                used.add(int(t))
                if t not in known:
                    problems.append(f"panel {num:02d}: tag E{t} is not an Entry in the notebook")
            if it["kind"] == "fig":
                f = ROOT / it["path"]
                if not f.exists():
                    problems.append(f"panel {num:02d}: missing figure {it['path']}")
                    continue
                if side:
                    fy += P.image(f, fx, fy, iw) + s["para_gap"]
                    continue
                y += P.image(f, ix, y, iw) + s["fig_gap"]
                if it["caption"]:
                    y = P.flow(it["caption"], ix, y, iw, s["caption_pt"], s["ink2"])
                y += s["para_gap"]
            elif it["kind"] == "bullet":
                ax.text(ix, y + s["body_pt"] * PT_MM, "•", fontsize=s["body_pt"], color=s["ink"],
                        family=s["font"], va="baseline", zorder=3)
                y = P.flow(it["text"], ix, y, iw, s["body_pt"], indent=8) + s["para_gap"]
            elif it["kind"] == "callout":
                y0c = y
                y = P.flow(it["text"], ix + 8, y + 7, iw - 16, s["callout_pt"]) + 7
                ax.add_patch(FancyBboxPatch((ix, y0c), iw, y - y0c, boxstyle="round,pad=0,rounding_size=3",
                                            fc=s["callout"], ec="none", zorder=2))
                y += s["para_gap"]
            else:
                y = P.flow(it["text"], ix, y, iw, s["body_pt"]) + s["para_gap"]
        if side:
            y = max(y, fy)
        y += s["pad"] - s["para_gap"]
        ax.add_patch(plt.Rectangle((x, card_top), w, y - card_top, fc=s["card"], ec="#dcdcdc", lw=1.2, zorder=1))
        y += s["panel_gap"]
        if p["wide"]:
            ys = [y] * s["cols"]
        else:
            count[col] += 1
            ys[col] = y + extra[col]
        if y - s["panel_gap"] > bottom + 0.5:
            problems.append(f"panel {num:02d} ({p['title'][:30]}) runs {y - s['panel_gap'] - bottom:.0f} mm "
                            f"into the footer")
    for k in ("body_pt", "caption_pt", "head_pt"):
        if s[k] < s["min_pt"]:
            problems.append(f"{k} = {s[k]} pt < min {s['min_pt']} pt")
    if not ASSETS.exists():
        print("note: poster/assets/ missing — logos skipped")

    left = [bottom - (v - s["panel_gap"] - extra[i]) for i, v in enumerate(ys)]
    if not save:
        plt.close(P.fig)
        return [lf / (n - 1) if n > 1 and lf > 0 else 0.0 for lf, n in zip(left, count)]
    OUT.mkdir(exist_ok=True)
    P.fig.savefig(OUT / "poster.pdf")
    P.fig.savefig(OUT / "poster.png", dpi=preview_dpi)
    print("column space distributed between panels (mm per gap): " +
          "  ".join(f"{i + 1}: {e:.0f}" for i, e in enumerate(extra)))
    print("entries cited: " + ", ".join(f"E{e}" for e in sorted(used)))
    print(f"wrote {OUT / 'poster.pdf'} and poster.png")
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for m in problems:
            print("  ! " + m)
        return 1
    print("checks passed")
    return 0


def text_height(P, items, ix, iw):
    """Height of the text items of a panel, measured without drawing (for vertical centring)."""
    s, y = SPEC, 0.0
    for it in items:
        if it["kind"] == "bullet":
            y = P.flow(it["text"], ix, y, iw, s["body_pt"], indent=8, draw=False) + s["para_gap"]
        elif it["kind"] == "callout":
            y = P.flow(it["text"], ix + 8, y + 7, iw - 16, s["callout_pt"], draw=False) + 7 + s["para_gap"]
        elif it["kind"] == "para":
            y = P.flow(it["text"], ix, y, iw, s["body_pt"], draw=False) + s["para_gap"]
    return y


if __name__ == "__main__":
    dpi = int(sys.argv[1]) if len(sys.argv) > 1 else 60      # optional preview dpi
    sys.exit(build(dpi, extra=build(save=False)))            # pass 1 measures, pass 2 balances the columns
