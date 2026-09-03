#!/usr/bin/env python3
"""Three explanatory flow charts for the thermodynamic-source method.

    fig1_reaction_selection.png   how a ModelSEED reaction acquires sources, and
                                  how TECRDB measurements are matched onto it
    fig2_grading.png              the gold/silver/bronze cascade
    fig3_recommendation.png       the per-target source-selection rule

These are DIAGRAMS, not data plots -- there is no data-to-position encoding. But
every count printed inside a box is READ FROM the shipped result files rather
than typed in, so the boxes cannot drift from the pipeline. The numbers are also
written to fig_flowchart_values.tsv beside the images, so prose that quotes a
box can be rechecked against a file.

Inputs (all tracked in git, all under reports/thermoSourceMethod/tables/):
    direction_maps_summary.json   per-source raw coverage, non-EMPTY total
    grade_calibration.json        vetoes, thresholds, anchor/proxy counts
    grade_frontier.tsv            grade counts per source
    recommendation_models.json    tolerances, per-target source mix, tau scales
Provenance of those files: ModelSEED dev @ 49563c6f (devsnap2), run 2026-08-12.
Regenerating them is NOT part of this script -- see reports/thermoSourceMethod/
THERMO_SOURCE_METHOD.md Part 12.

DESIGN NOTES (the things a change request will touch)
  * Node role is encoded by BORDER hue on a near-white fill, not by a saturated
    fill. Flow-chart nodes are mostly text, and text needs a light background;
    a saturated fill would force either white text (fails on the yellow) or a
    contrast warning.
  * Role palette validated with the dataviz validator in the legend order used
    below: worst adjacent CVD dE 15.3 (deutan), normal-vision 20.8. Reject nodes
    additionally carry a DASHED border, so role never rests on hue alone.
  * Layout is hand-placed on a 0-100 grid per figure. There is no auto-layout;
    if you add a node, move its neighbours.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ANALYSIS_DIR = Path(os.environ.get("CORE_MODELS_ANALYSIS_DIR",
                                   "/scratch/ctaylor/core_models_analysis"))
TABLES = Path(os.environ.get("FLOWCHART_TABLES",
                             str(ANALYSIS_DIR / "reports" / "thermoSourceMethod" / "tables")))
OUT_DIR = Path(os.environ.get("FLOWCHART_OUT",
                              str(ANALYSIS_DIR / "reports" / "thermoSourceMethod" / "figures")))

# --- ink + surface tokens, shared with plot_graded_fba.py so the report reads as one set
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#898781"
SURFACE, RULE = "#fcfcfb", "#c3c2b7"
# --- node roles. Legend order is the validated order; do not reshuffle without re-running
#     dataviz/scripts/validate_palette.js.
ROLE = {
    "data":     "#2a78d6",   # something read from disk
    "outcome":  "#1baf7a",   # a label or choice the pipeline emits
    "process":  "#4a3aa7",   # a computation
    "decision": "#eda100",   # a branch
    "reject":   "#e34948",   # dropped / abstained -- also drawn DASHED
}
TINT = 0.90          # how far each border hue is mixed toward white for the fill
FS_TITLE, FS_NODE, FS_SUB, FS_EDGE = 13.0, 9.3, 8.2, 8.0


def _tint(hex_colour: str, amount: float = TINT) -> tuple:
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return tuple(c + (1.0 - c) * amount for c in (r, g, b))


class Chart:
    """A 0-100 x 0-100 canvas with boxes and arrows."""

    def __init__(self, w, h, title, subtitle=""):
        self.fig, self.ax = plt.subplots(figsize=(w, h))
        self.fig.patch.set_facecolor(SURFACE)
        self.ax.set_facecolor(SURFACE)
        self.ax.set_xlim(0, 100)
        self.ax.set_ylim(0, 100)
        self.ax.axis("off")
        self.ax.set_title(title, color=INK, fontsize=FS_TITLE, loc="left", pad=34,
                          fontweight="semibold")
        if subtitle:
            self.ax.text(0, 101.6, subtitle, color=INK2, fontsize=FS_SUB + 0.6,
                         va="bottom", ha="left", linespacing=1.45)
        self.boxes = {}

    def box(self, key, x, y, w, h, label, sub="", role="process"):
        """x, y = CENTRE. Returns nothing; look the node up by key for arrows.

        The box AUTO-GROWS if the text will not fit in the requested height.
        Hand-tuning every height is how flow charts rot: someone adds a line to
        a label months later and the text silently spills over the border.
        """
        # Lay label and sub out as one vertically centred block. Doing this by
        # eye breaks the moment a label wraps to two lines, which is most of them.
        lh_l, lh_s, gap, pad = 2.85, 2.45, 0.55, 1.8
        n_l = label.count("\n") + 1
        n_s = (sub.count("\n") + 1) if sub else 0
        total = n_l * lh_l + (gap + n_s * lh_s if n_s else 0.0)
        h = max(h, total + pad)
        colour = ROLE[role]
        p = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                           boxstyle="round,pad=0,rounding_size=1.6",
                           linewidth=1.7, edgecolor=colour, facecolor=_tint(colour),
                           linestyle=(0, (4, 2.2)) if role == "reject" else "solid",
                           zorder=3)
        self.ax.add_patch(p)
        top = y + total / 2.0
        self.ax.text(x, top - n_l * lh_l / 2.0, label, ha="center", va="center",
                     fontsize=FS_NODE, color=INK, zorder=4, linespacing=1.42)
        if sub:
            self.ax.text(x, top - n_l * lh_l - gap - n_s * lh_s / 2.0, sub,
                         ha="center", va="center", fontsize=FS_SUB, color=INK2,
                         zorder=4, linespacing=1.38)
        self.boxes[key] = (x, y, w, h)

    def _port(self, key, side):
        x, y, w, h = self.boxes[key]
        return {"t": (x, y + h / 2), "b": (x, y - h / 2),
                "l": (x - w / 2, y), "r": (x + w / 2, y)}[side]

    def arrow(self, a, b, sa="b", sb="t", label="", rad=0.0, lx=0, ly=0, dim=False):
        p0, p1 = self._port(a, sa), self._port(b, sb)
        self.ax.add_patch(FancyArrowPatch(
            p0, p1, arrowstyle="-|>", mutation_scale=11,
            linewidth=1.25, color=RULE if dim else INK3,
            connectionstyle=f"arc3,rad={rad}", zorder=2,
            shrinkA=1.5, shrinkB=2.5))
        if label:
            self.ax.text((p0[0] + p1[0]) / 2 + lx, (p0[1] + p1[1]) / 2 + ly, label,
                         ha="center", va="center", fontsize=FS_EDGE, color=INK2,
                         zorder=5, bbox=dict(boxstyle="round,pad=0.22", fc=SURFACE,
                                             ec="none"))

    def legend(self, roles, y=-4.5):
        x = 0.0
        for role, text in roles:
            c = ROLE[role]
            self.ax.add_patch(FancyBboxPatch(
                (x, y - 1.0), 2.4, 2.0, boxstyle="round,pad=0,rounding_size=0.7",
                linewidth=1.6, edgecolor=c, facecolor=_tint(c),
                linestyle=(0, (3, 1.8)) if role == "reject" else "solid",
                clip_on=False, zorder=3))
            self.ax.text(x + 3.4, y, text, fontsize=FS_SUB, color=INK2,
                         va="center", ha="left", clip_on=False)
            x += 4.2 + len(text) * 1.02

    def note(self, text, y=-9.5):
        self.ax.text(0, y, text, fontsize=FS_SUB, color=INK3, va="top", ha="left",
                     clip_on=False, linespacing=1.5)

    def save(self, name, bottom=0.16):
        self.fig.subplots_adjust(left=0.015, right=0.985, top=0.88, bottom=bottom)
        self.fig.savefig(OUT_DIR / name, dpi=200, facecolor=SURFACE)
        plt.close(self.fig)
        print(f"  wrote {name}")


# ------------------------------------------------------------------ the inputs
def load_values() -> dict:
    dm = json.load(open(TABLES / "direction_maps_summary.json"))
    cal = json.load(open(TABLES / "grade_calibration.json"))
    rec = json.load(open(TABLES / "recommendation_models.json"))
    fr = {r["source"]: r for r in csv.DictReader(
        open(TABLES / "grade_frontier.tsv"), delimiter="\t")}
    v = {
        "n_nonempty": dm["n_reactions_nonempty"],
        "cov_gc": dm["coverage"]["gc"], "cov_eq": dm["coverage"]["eq"],
        "cov_dg": dm["coverage"]["dgpms"], "n_feasible": dm["coverage"]["graded"],
        "veto_sent": cal["vetoes"]["eq_sentinel"],
        "veto_mnx": cal["vetoes"]["eq_mnx_collision"],
        "veto_quin": cal["vetoes"]["dgpms_quinone"],
        "p_gold": cal["thresholds"]["p_gold"], "p_silver": cal["thresholds"]["p_silver"],
        "r_corrob": cal["thresholds"]["r_corrob"], "z_corrob": cal["thresholds"]["z_corrob"],
        "r_out": cal["thresholds"]["r_outvote"], "z_out": cal["thresholds"]["z_outvote"],
        "m_gold": cal["thresholds"]["meas_gold"], "m_silver": cal["thresholds"]["meas_silver"],
        "anchor": cal["p_ok_models"]["GC"]["n_anchor"],
        "proxy_gc": cal["p_ok_models"]["GC"]["n_proxy"],
        "proxy_eq": cal["p_ok_models"]["EQ"]["n_proxy"],
        "proxy_dg": cal["p_ok_models"]["DGPMS"]["n_proxy"],
        "tol_dir": rec["targets"]["direction"]["tolerance"],
        "tol_mag": rec["targets"]["magnitude"]["tolerance"],
        "kept_dir": rec["targets"]["direction"]["n_kept"],
        "kept_mag": rec["targets"]["magnitude"]["n_kept"],
    }
    for src, key in [("Group contribution", "gc"), ("eQuilibrator", "eq"),
                     ("dGPredictor-ModelSEED", "dg"), ("TECRDB", "tec")]:
        for tier in ("gold", "silver", "bronze"):
            v[f"{tier}_{key}"] = int(fr[src][f"n_{tier}"])
        v[f"graded_{key}"] = int(fr[src]["n_graded"])
    mix = rec["targets"]["direction"]["mix"]
    for lab, key in [("eQuilibrator", "eq"), ("dGPredictor-ModelSEED", "dg"),
                     ("TECRDB", "tec"), ("Group contribution", "gc")]:
        v[f"dirmix_{key}"] = mix.get(lab, 0)
    v["skeleton"] = v["silver_tec"]
    v["no_source"] = v["n_nonempty"] - v["n_feasible"]
    v["abstain_dir"] = v["n_feasible"] - v["kept_dir"]
    v["abstain_mag"] = v["n_feasible"] - v["kept_mag"]
    return v


def n(x) -> str:
    return f"{x:,}"


# -------------------------------------------------------------------- figure 1
def fig_selection(v):
    c = Chart(13.6, 8.0,
              "1 · Which reactions get thermodynamic data, and how experiments are matched onto them",
              "Every ModelSEED reaction enters at the top left. TECRDB measurements enter at the right "
              "and are matched by structure, never by identifier.")
    c.box("all", 25, 92, 32, 8, "ModelSEED  dev @ 49563c6f",
          f"{n(v['n_nonempty'])} non-EMPTY reactions", "data")
    c.box("gc", 9, 73, 16, 12, "Group\nContribution", n(v["cov_gc"]), "data")
    c.box("eq", 25.5, 73, 16, 12, "eQuilibrator", n(v["cov_eq"]), "data")
    c.box("dg", 42, 73, 16, 12, "dGPredictor\n-ModelSEED", n(v["cov_dg"]), "data")
    for k in ("gc", "eq", "dg"):
        c.arrow("all", k, "b", "t", rad=0.0 if k == "gc" else 0.03)

    c.box("veto", 25.5, 55, 49, 11,
          "Drop values that are not estimates",
          f"eQuilibrator sentinel σ>100  {n(v['veto_sent'])}      "
          f"MetaNetX collision  {n(v['veto_mnx'])}      "
          f"dGPredictor on quinones  {n(v['veto_quin'])}", "reject")
    for k in ("gc", "eq", "dg"):
        c.arrow(k, "veto", "b", "t", rad=0.0)

    c.box("feas", 12, 36.5, 21, 12, "Feasible sources",
          f"{n(v['n_feasible'])} reactions", "outcome")
    c.box("none", 38, 36.5, 24, 12, "No source at all",
          f"{n(v['no_source'])} reactions\na hard ceiling", "reject")
    c.arrow("veto", "feas", "b", "t", rad=-0.12, label="≥1 survives", ly=1.2)
    c.arrow("veto", "none", "b", "t", rad=0.12, label="none", ly=1.2)

    # --- TECRDB branch
    c.box("tec", 78, 92, 42, 8, "TECRDB  (Zenodo 3978440)",
          "4,544 measured K′ rows", "data")
    c.box("struct", 78, 73.5, 42, 8.5, "Resolve every compound to an RDKit structure key",
          "", "process")
    c.box("match", 78, 55, 42, 13,
          "Match as a reaction, not an ID",
          "(reactant multiset, product multiset)\nprotons dropped · both directions", "process")
    c.arrow("tec", "struct", "b", "t")
    c.arrow("struct", "match", "b", "t")

    c.box("stereo", 67, 34.5, 21, 13, "stereo_exact",
          f"full InChIKey · {n(v['anchor'])}\ndistinguishes anomers, D/L", "outcome")
    c.box("skel", 89.5, 34.5, 19, 13, "skeleton",
          f"connectivity only · {n(v['skeleton'])}\ncan conflate stereoisomers", "decision")
    c.arrow("match", "stereo", "b", "t", rad=-0.1)
    c.arrow("match", "skel", "b", "t", rad=0.1)

    c.box("anchor", 67, 12, 21, 13, "THE ANCHOR SET",
          f"{n(v['anchor'])} reactions where\nthe truth is known", "outcome")
    c.arrow("stereo", "anchor", "b", "t")
    c.box("cap", 89.5, 12, 19, 13, "Used, capped\nat SILVER",
          "measurement sound,\nmatch may not be", "decision")
    c.arrow("skel", "cap", "b", "t")

    c.legend([("data", "read from disk"), ("process", "computation"),
              ("decision", "branch"), ("outcome", "pipeline output"),
              ("reject", "dropped")], y=-3.0)
    c.note("The 802 anchor reactions are a subset of the 33,289 feasible ones — the two columns meet "
           "there. That anchor set is what every\ncalibration and every validation in figures 2 and 3 "
           "is fitted and scored on. It is 1.4% of the database, and all of it is\nwell-studied "
           "central metabolism: the easy part.", y=-7.0)
    c.save("fig4_reaction_selection.png", bottom=0.19)


# -------------------------------------------------------------------- figure 2
def fig_grading(v):
    c = Chart(13.6, 10.6, "2 · The grading algorithm — a trust label for every source, on every reaction",
              "Applied independently to each source, so on one reaction eQuilibrator can be GOLD "
              "while Group Contribution is BRONZE.")
    # `in` is centred over BOTH downstream columns so the two arrows leave from
    # one point and diverge; routing one of them past UNGRADED made it read as
    # though UNGRADED fed STEP 2.
    c.box("in", 44, 93, 46, 9, "One (reaction, source) pair",
          "its ΔG′° and its own reported σ", "data")
    c.box("ung", 86, 93, 26, 9, "UNGRADED",
          "no value, or vetoed in figure 1", "reject")
    c.arrow("in", "ung", "r", "l", label="not\nfeasible", ly=3.6)

    c.box("cal", 22, 76, 38, 13, "STEP 1 · Calibrate σ",
          f"isotonic fit, per source, of σ against real error\n"
          f"anchor {n(v['anchor'])} at weight 3  +  proxy "
          f"{n(v['proxy_eq'])}–{n(v['proxy_dg'])} at weight 1", "process")
    c.box("fuse", 72, 76, 42, 13, "STEP 2 · Compare the sources",
          "weight each by ê, then ask whether the spread\n"
          "is what their σ predicts → Birge R, residual z", "process")
    c.arrow("in", "cal", "b", "t", rad=0.10)
    c.arrow("in", "fuse", "b", "t", rad=-0.10)

    c.box("out1", 22, 62, 38, 7.5, "ê  expected error       p  P(error ≤ 2)", "", "outcome")
    c.box("out2", 72, 62, 42, 7.5, "R  do they agree       z  who is the outlier", "", "outcome")
    c.arrow("cal", "out1", "b", "t")
    c.arrow("fuse", "out2", "b", "t")

    y0, dy, bw, bx = 52, 11.3, 84, 47
    c.box("r2", bx, y0, bw, 8.6,
          f"RULE 2 · own confidence          p ≥ {v['p_gold']} → GOLD          "
          f"p ≥ {v['p_silver']} → SILVER          else BRONZE", "", "decision")
    c.box("r3", bx, y0 - dy, bw, 8.6,
          f"RULE 3 · corroborated             BRONZE and R ≤ {v['r_corrob']} and "
          f"z ≤ {v['z_corrob']:.0f}  →  SILVER,  never higher", "", "decision")
    c.box("r4", bx, y0 - 2 * dy, bw, 8.6,
          f"RULE 4 · outvoted                    R > {v['r_out']:.0f} and "
          f"z > {v['z_out']:.0f}  →  one tier down", "", "decision")
    c.box("r1", bx, y0 - 3 * dy, bw, 10.5,
          f"RULE 1 · measured                  |error| ≤ {v['m_gold']:.0f} → GOLD      "
          f"≤ {v['m_silver']:.0f} → SILVER      else BRONZE",
          "applied last, so it overrides everything above", "outcome")
    c.arrow("out1", "r2", "b", "t", rad=-0.06)
    c.arrow("out2", "r3", "b", "t", rad=0.16)
    c.arrow("r2", "r3", "b", "t")
    c.arrow("r3", "r4", "b", "t")
    c.arrow("r4", "r1", "b", "t")

    gy, gw = 6.5, 21
    c.box("g_tec", 12, gy, gw, 9, "TECRDB",
          f"{n(v['gold_tec'])} gold · {n(v['silver_tec'])} silver", "outcome")
    c.box("g_eq", 35.5, gy, gw, 9, "eQuilibrator",
          f"{n(v['gold_eq'])} / {n(v['silver_eq'])} / {n(v['bronze_eq'])}", "outcome")
    c.box("g_dg", 59, gy, gw, 9, "dGPredictor-MS",
          f"{n(v['gold_dg'])} / {n(v['silver_dg'])} / {n(v['bronze_dg'])}", "outcome")
    c.box("g_gc", 82.5, gy, gw, 9, "Group Contribution",
          f"{n(v['gold_gc'])} / {n(v['silver_gc'])} / {n(v['bronze_gc'])}", "outcome")
    c.arrow("r1", "g_eq", "b", "t", rad=0.0, dim=True)

    c.legend([("data", "input"), ("process", "computation"),
              ("decision", "cascade rule"), ("outcome", "output"),
              ("reject", "dropped")], y=-5.0)
    c.note("Bottom row is gold / silver / bronze per source, database-wide. TECRDB is graded "
           "separately — gold everywhere it exists, except\nskeleton-tier matches, capped at silver. "
           "Group Contribution reaches gold only through RULE 1: its own confidence never\nclears "
           "0.90 anywhere in the database, so all 309 of its golds are measured reactions.", y=-9.5)
    c.save("fig5_grading.png", bottom=0.21)


# -------------------------------------------------------------------- figure 3
def fig_recommendation(v):
    c = Chart(13.6, 8.6, "3 · The recommendation algorithm — which single source to actually use",
              "A grade says how much to trust a number; a recommendation picks one. "
              "The rule depends on what the number is for.")
    c.box("in", 20, 93, 28, 7.5, "One reaction",
          f"{n(v['n_feasible'])} have at least one feasible source", "data")
    c.box("meas", 20, 78, 28, 8.5, "Is it measured?",
          "TECRDB present → return it", "decision")
    c.arrow("in", "meas", "b", "t")
    c.box("tec_out", 62, 78, 22, 8.5, "Use the experiment",
          f"{n(v['dirmix_tec'])} reactions", "outcome")
    c.arrow("meas", "tec_out", "r", "l", label="yes", ly=1.6)

    c.box("feas", 20, 62, 28, 8.5, "Drop vetoed sources",
          "sentinel · collision · quinone", "reject")
    c.arrow("meas", "feas", "b", "t", label="no", lx=-4.2)

    c.box("tau", 20, 46.5, 28, 8.5, "Calibrate σ → τ",
          f"τ = k·ê/√(2/π),  k = {v['tau_kgc'] if 'tau_kgc' in v else '0.54–0.65'}", "process")
    c.arrow("feas", "tau", "b", "t")

    c.box("split", 20, 31, 28, 8.5, "What is the number FOR?", "", "decision")
    c.arrow("tau", "split", "b", "t")

    c.box("mag", 55, 42, 30, 11, "MAGNITUDE  →  use the uncertainty",
          f"pick argmin ê\nabstain if ê > {v['tol_mag']:.0f} kcal/mol", "process")
    c.box("dir", 55, 20, 30, 11, "DIRECTION  →  do NOT use it to choose",
          f"fixed priority  eQ ▸ dGPredictor ▸ GC\nabstain if risk > {v['tol_dir']}", "process")
    c.arrow("split", "mag", "r", "l", rad=-0.18, label="the energy", ly=2.2)
    c.arrow("split", "dir", "r", "l", rad=0.10, label="the direction", ly=-2.2)

    c.box("mag_out", 89, 42, 20, 11, f"{n(v['kept_mag'])} answered",
          f"{n(v['abstain_mag'])} abstained\n(64% of the covered set)", "outcome")
    c.box("dir_out", 89, 20, 20, 11, f"{n(v['kept_dir'])} answered",
          f"{n(v['abstain_dir'])} abstained\n(18% of the covered set)", "outcome")
    c.arrow("mag", "mag_out", "r", "l")
    c.arrow("dir", "dir_out", "r", "l")

    c.box("why", 20, 12, 30, 17, "Why not use σ to choose\nthe direction source?",
          "every uncertainty-based rule lost to\nthe flat priority list:  argmin risk 90.9%,\n"
          "argmin ê 93.7%,  fixed priority 95.9%", "reject")
    c.arrow("dir", "why", "l", "r", rad=0.0, dim=True)

    c.legend([("data", "input"), ("process", "computation"),
              ("decision", "branch"), ("outcome", "output"),
              ("reject", "dropped / cautionary")], y=-4.0)
    c.note("The risk is P(this source's own call survives its own uncertainty) — precision, not "
           "accuracy. A confidently wrong\nsource scores well on it. So it is used to decide whether "
           "to answer at all, never to decide which source answers.", y=-8.0)
    c.save("fig6_recommendation.png", bottom=0.20)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    v = load_values()
    fig_selection(v)
    fig_grading(v)
    fig_recommendation(v)
    stats = OUT_DIR / "fig_flowchart_values.tsv"
    with open(stats, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["key", "value", "note"])
        w.writerow(["_source", str(TABLES), "every value below is read from these tables"])
        for k in sorted(v):
            w.writerow([k, v[k], ""])
    print(f"  wrote {stats.name} ({len(v)} values)")


if __name__ == "__main__":
    main()
