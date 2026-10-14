r"""poster/poster.md -> poster/out/poster.pdf (+ PNG preview). Run: python poster/build.py [preview_dpi]

Same idea as a print pipeline: the manuscript is canonical and read-only, every coordinate comes
from SPEC (mm, page origin top-left), and the build ends with machine checks. Looking at the PNG
is still required; the checks only catch what a machine can see.

Look (v3, lab template). Matches the CUPT lab posters for KPS 2026 (J. Jang, poster_refs): navy
header with the IHS shield left and the white CUPT logo right, three-colour stripe, white page,
thin dark-outlined boxes with "■ Title" heads, an Acknowledgement box at the bottom.
Text in Source Sans 3 (the lab posters use Arial). A QR code to the repository sits in the
Acknowledgement box (poster/assets/qr_repo.png, made once with segno).
v2 (gray page, shadowed cards, orange number badges, Plex Sans) is in git history (c381842).
Kept from v2:
- text in IBM Plex Sans (poster/fonts, SIL OFL); math in STIX, so a variable is a serif italic
  and is told apart from the sans body at a glance. Plex Sans was the math face too until the
  specimen showed what it costs: it has no \cdot and no \sim, so "h . u_0" and "R ~ 1" were
  silently falling back to DejaVu, and its \int is a text-sized glyph that vanishes in the
  Theta_n equation. Computer Modern ("cm") fixes both but its display integral overshoots the
  line box and collides with the line below, so STIX is the one that fits this layout. Set
  TRATT_MATH=cm|plex|dejavuserif to re-compare. Absolute values are \left|..\right|: a bare
  STIX | is the height of a lowercase l and reads as one;
- bullets are small squares at x-height; the key message is a beige band (as on the ICAMD poster);
- "%%" stat tiles for the few numbers a passer-by should take away.

Manuscript syntax
  front matter: title ("//" = line break), authors, affiliation, contact, footer
  "## Title"      a panel (numbered automatically)        "@wide"   next panel spans both columns
  "@column"       continue in the next column             "> text"  key-message band
  "- text"        bullet                                  "![caption](path)" figure (width = card)
  "%% value | label"  stat tile (consecutive lines form one row)
  {E12}           notebook entry tag (stripped from print, required on every item)

Checks
1. every paragraph, bullet, band, tile and figure carries an {E..} tag, and each tagged Entry exists
   in notes/lab_notebook.md (nothing goes on the poster without its entry);
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
from matplotlib import font_manager  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
for _f in sorted((ROOT / "poster" / "fonts").glob("*.ttf")):
    font_manager.fontManager.addfont(str(_f))
FONT = "Source Sans 3"   # poster/fonts, SIL OFL; humanist, narrower than Arial, reads well at distance
# Word-by-word placement needs text widths that scale linearly with dpi; hinting breaks that and
# makes long words overrun the next one in a preview rendered at another dpi.
import os  # noqa: E402

# Math face. IBM Plex Sans has no \cdot and no \sim, so a custom all-Plex fontset silently
# falls back to DejaVu for "h . u_0" and "R ~ 1"; its \int is also a text-sized glyph, which
# the Theta_n equation needs to be bigger than. Set TRATT_MATH to compare.
MATH = os.environ.get("TRATT_MATH", "stix")
_math = ({"mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic",
          "mathtext.bf": f"{FONT}:bold", "mathtext.sf": FONT}
         if MATH == "plex" else {"mathtext.fontset": MATH})
matplotlib.rcParams.update({"text.hinting": "none", "text.hinting_factor": 1, "font.family": FONT,
                            "pdf.fonttype": 42,   # TrueType in the PDF: editable, accepted by print shops
                            **_math})
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

MS = ROOT / "poster" / "poster.md"
OUT = ROOT / "poster" / "out"
ASSETS = ROOT / "poster" / "assets"
NOTEBOOK = ROOT / "notes" / "lab_notebook.md"

# A0 portrait, same as the ICAMD 2025 poster (841 x 1189 mm). Check the KPS rule before printing.
SPEC = dict(
    W=841, H=1189, margin=16, gap=14, cols=2, pad=13,
    header_h=130, stripe_h=5, footer_h=56,
    title_pt=52, authors_pt=30, head_pt=40, badge_pt=30, body_pt=26, callout_pt=29, caption_pt=21,
    stat_pt=50, stat_label_pt=21, min_pt=20,
    leading=1.32, para_gap=7, panel_gap=11, fig_gap=4,
    navy="#0f2a3a", accent="#e8743b", page="#ffffff", card="#ffffff", edge="#4a4a4a", key="#fbf1d9",
    ink="#111111", ink2="#4f5560", tint="#f2f4f7",
    stripe=("#3a3a3a", "#e04b3a", "#f08a3c"),
)
PT_MM = 25.4 / 72
CAP, XH = 0.66, 0.49        # Source Sans 3 cap height and x-height / em


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
        if block.startswith("%%"):
            tiles = []
            for line in block.splitlines():
                v, lab = line[2:].split("|", 1)
                tiles.append(dict(value=v.strip(), label=TAG.sub("", lab).strip(), tags=TAG.findall(lab)))
            cur["items"].append(dict(kind="stats", tiles=tiles, tags=[t for x in tiles for t in x["tags"]],
                                     ok=all(x["tags"] for x in tiles)))
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

    def width(self, tok, pt, weight="normal"):
        key = (tok, pt, weight)
        if key not in self._w:
            t = self.ax.text(0, 0, tok, fontsize=pt, fontweight=weight, family=FONT)
            self._w[key] = t.get_window_extent(self.r).width / self.fig.dpi * 25.4
            t.remove()
        return self._w[key]

    def flow(self, text, x, y, w, pt, color=None, weight="normal", indent=0.0, align="left", z=3, draw=True,
             strong=None):
        """Word-by-word placement with wrapping; $..$ and **..** stay atomic, **..** is semibold (in colour
        `strong` if given). Returns the y (mm) just below the last line."""
        color = color or SPEC["ink"]
        lh = pt * PT_MM * SPEC["leading"]
        sp = self.width("a a", pt) - 2 * self.width("a", pt)
        y = y + pt * PT_MM                                  # baseline of the first line
        toks, prev = [], 0
        for m in TOKEN.finditer(text):
            glue = m.start() == prev and bool(toks)
            prev = m.end()
            s, wgt, col = m[0], weight, color
            if s.startswith("**"):
                s, wgt, col = s[2:-2], "semibold", strong or color
            toks.append((s, wgt, col, glue))
        lines, line, lw = [], [], 0.0
        for s, wgt, col, g in toks:
            tw = self.width(s, pt, wgt)
            add = tw + (0 if g or not line else sp)
            if line and lw + add > w - indent:
                k = len(line)                                # a glued run ending the line moves along
                while g and k > 1 and line[k - 1][3]:
                    k -= 1
                if g:
                    k -= 1
                carry = line[k:] if k > 0 else []
                kept = line[:k] if k > 0 else line
                lines.append((kept, sum(t[4] for t in kept) + sp * sum(1 for t in kept[1:] if not t[3])))
                line = [(cs, cw, cc, i > 0 and cg, ctw) for i, (cs, cw, cc, cg, ctw) in enumerate(carry)]
                lw = sum(t[4] for t in line) + sp * sum(1 for t in line[1:] if not t[3])
                add = tw if not line else tw + (0 if g else sp)
            line.append((s, wgt, col, g if line else False, tw))
            lw += add
        if line:
            lines.append((line, lw))
        for line, lw in lines:
            cx = x + indent + ((w - indent - lw) / 2 if align == "center" else 0)
            first = True
            for s, wgt, col, g, tw in line:
                if not g and not first:
                    cx += sp
                if draw:
                    self.ax.text(cx, y, s, fontsize=pt, color=col, family=FONT, zorder=z, fontweight=wgt,
                                 va="baseline")
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

    def box(self, x, y, w, h, fc, r=2.5, z=1, ec="none", lw=0.0):
        self.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                         fc=fc, ec=ec, lw=lw, zorder=z))


# ───────────── items ─────────────

def bullet(P, it, x, y, w, draw=True):
    s = SPEC
    if draw:
        d = s["body_pt"] * PT_MM * 0.24                     # marker side
        mid = y + s["body_pt"] * PT_MM * (1 - XH / 2)       # middle of the x-height of the first line
        P.box(x + 0.5, mid - d / 2, d, d, s["ink"], r=0, z=3)
    return P.flow(it["text"], x, y, w, s["body_pt"], indent=7, draw=draw,
                  strong=s["accent"]) + s["para_gap"]


def callout(P, it, x, y, w, draw=True):
    s = SPEC
    y1 = P.flow(it["text"], x + 10, y + 9, w - 20, s["callout_pt"], s["ink"], draw=draw, strong=s["accent"]) + 8
    if draw:
        P.box(x, y, w, y1 - y, s["key"], r=0, z=2)
    return y1 + s["para_gap"]


def stats(P, it, x, y, w, draw=True):
    s = SPEC
    n = len(it["tiles"])
    g = 6
    tw = (w - (n - 1) * g) / n
    vh = s["stat_pt"] * PT_MM * 1.05
    bottoms = []
    for i, t in enumerate(it["tiles"]):
        tx = x + i * (tw + g)
        if draw:
            P.flow(t["value"], tx, y + 6, tw, s["stat_pt"], s["accent"], weight="semibold", align="center")
        bottoms.append(P.flow(t["label"], tx + 4, y + 6 + vh + 3, tw - 8, s["stat_label_pt"], s["ink2"],
                              align="center", draw=draw))
    y1 = max(bottoms) + 6
    if draw:
        for i in range(n):
            P.box(x + i * (tw + g), y, tw, y1 - y, s["tint"], r=0, z=2)
    return y1 + s["para_gap"]


def text_height(P, items, ix, iw):
    """Height of the text items of a panel, measured without drawing (for vertical centring)."""
    s, y = SPEC, 0.0
    for it in items:
        if it["kind"] == "bullet":
            y = bullet(P, it, ix, y, iw, draw=False)
        elif it["kind"] == "callout":
            y = callout(P, it, ix, y, iw, draw=False)
        elif it["kind"] == "stats":
            y = stats(P, it, ix, y, iw, draw=False)
        elif it["kind"] == "para":
            y = P.flow(it["text"], ix, y, iw, s["body_pt"], draw=False,
                      strong=s["accent"]) + s["para_gap"]
    return y


# ───────────── page ─────────────

def build(preview_dpi=60, extra=None, save=True):
    """One layout pass. extra[c] is added after each panel of column c except the last, so the
    columns end level with the footer (filled in by a first, unsaved pass; see __main__)."""
    s = SPEC
    extra = extra or [0.0] * s["cols"]
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
        logo_w = 100
        P.image(ASSETS / "logo_ihs.png", M + 6, 22, 80 * 0.7, h=80 * 1.0)
    if (ASSETS / "logo_cupt_white.png").exists():
        P.image(ASSETS / "logo_cupt_white.png", W - M - 100, 40, 100)
    tw = W - 2 * M - 2 * logo_w
    y = 20 - s["title_pt"] * PT_MM * (s["leading"] - 1)
    for line in meta.get("title", "").split("//"):           # "//" in the title = manual line break
        y = P.flow(line.strip(), M + logo_w, y + s["title_pt"] * PT_MM * (s["leading"] - 1), tw, s["title_pt"],
                   "white", weight="bold", align="center") - s["title_pt"] * PT_MM * 0.35
    y = P.flow(meta.get("authors", ""), M + logo_w, y + 13, tw, s["authors_pt"], "white", align="center")
    P.flow(meta.get("affiliation", ""), M + logo_w, y + 4, tw, s["authors_pt"] - 8, "#c9d6e0", align="center")

    # footer: an Acknowledgement box, as on the lab posters
    fy = H - s["footer_h"]
    ax.text(M + s["pad"], fy + 6 + s["head_pt"] * PT_MM * 0.9, "■ Acknowledgement", fontsize=s["head_pt"] - 6,
            color=s["ink"], family=FONT, fontweight="bold", va="baseline", zorder=3)
    qr = ASSETS / "qr_repo.png"
    qw = s["footer_h"] - M * 0.75 - 6 if qr.exists() else 0
    if qw:
        P.image(qr, W - M - 3 - qw, fy + 3, qw, h=qw)
        P.flow("code and lab notebook", W - M - 6 - qw - 70, fy + 3 + qw / 2 - 4, 68, s["caption_pt"],
               s["ink2"], align="center")
    P.flow(" · ".join(v for k, v in meta.items() if k in ("footer", "contact")), M + s["pad"],
           fy + 8 + s["head_pt"] * PT_MM, W - 2 * M - 2 * s["pad"] - (qw + 80 if qw else 0), s["caption_pt"] + 2,
           s["ink"])
    P.box(M, fy, W - 2 * M, s["footer_h"] - M * 0.75, s["card"], r=0, z=1, ec=s["edge"], lw=1.6)

    cw = (W - 2 * M - (s["cols"] - 1) * s["gap"]) / s["cols"]
    top = s["header_h"] + s["stripe_h"] + 2 + s["panel_gap"] + 2
    bottom = fy - 4
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
        card_top = y0
        # section head inside the box: "■ Title", as on the lab posters
        y = P.flow("■ " + p["title"], ix, y0 + s["pad"] * 0.7, iw, s["head_pt"], s["ink"], weight="bold") + 8
        # A wide panel with a figure is laid out side by side: text left, figure right (half width each).
        side = p["wide"] and any(it["kind"] == "fig" for it in p["items"])
        if side:
            fx, ix, iw = ix + iw / 2 + s["gap"] / 2, ix, iw / 2 - s["gap"] / 2
            fh = sum(iw * (lambda a: a.shape[0] / a.shape[1])(mpimg.imread(ROOT / it["path"]))
                     for it in p["items"] if it["kind"] == "fig" and (ROOT / it["path"]).exists())
            th = text_height(P, p["items"], ix, iw) - s["para_gap"]
            fyy = y + max(0.0, (th - fh) / 2)
            y += max(0.0, (fh - th) / 2)
        for it in p["items"]:
            if not it["tags"] or (it["kind"] == "stats" and not it["ok"]):
                problems.append(f"panel {num:02d}: untagged {it['kind']}: {it.get('text', it.get('caption', ''))[:50]!r}")
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
                    fyy += P.image(f, fx, fyy, iw) + s["para_gap"]
                    continue
                y += P.image(f, ix, y, iw) + s["fig_gap"]
                if it["caption"]:
                    y = P.flow(it["caption"], ix, y, iw, s["caption_pt"], s["ink2"],
                                strong=s["accent"])
                y += s["para_gap"]
            elif it["kind"] == "bullet":
                y = bullet(P, it, ix, y, iw)
            elif it["kind"] == "callout":
                y = callout(P, it, ix, y, iw)
            elif it["kind"] == "stats":
                y = stats(P, it, ix, y, iw)
            else:
                y = P.flow(it["text"], ix, y, iw, s["body_pt"],
                           strong=s["accent"]) + s["para_gap"]
        if side:
            y = max(y, fyy)
        y += s["pad"] - s["para_gap"]
        P.box(x, card_top, w, y - card_top, s["card"], r=0, z=1, ec=s["edge"], lw=1.6)
        y += s["panel_gap"]
        if p["wide"]:
            ys = [y] * s["cols"]
        else:
            count[col] += 1
            ys[col] = y + extra[col]
        if y - s["panel_gap"] > bottom + 0.5:
            problems.append(f"panel {num:02d} ({p['title'][:30]}) runs {y - s['panel_gap'] - bottom:.0f} mm "
                            f"into the footer")
    for k in ("body_pt", "caption_pt", "head_pt", "stat_label_pt"):
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


if __name__ == "__main__":
    dpi = int(sys.argv[1]) if len(sys.argv) > 1 else 60      # optional preview dpi
    sys.exit(build(dpi, extra=build(save=False)))            # pass 1 measures, pass 2 balances the columns
