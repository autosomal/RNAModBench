# Figure S5 (rebuilt 2026-09-21) — legend, paste-ready

Figure file: `FigureS5_rev.pdf` (**240 × 175 mm**, drawn 1:1 at
the printed size so the text stays large relative to the figure; Arial embedded,
no gridlines, no in-panel annotation text, **the 13 tool names slanted 45° as
the only rotated text**, bold panel letters A–E printed once per row).
Produced by `figures/figureS6/src/62_figS5_figure.py` from the frozen
tables of `40_fig6_combination.py` and `61_figS5_sitequality.py` (nothing is
recomputed in the figure script) and verified by `63_verify_figS5.py`
(ALL CHECKS PASSED).  Tables live in
`figures/figure6/tables/` (`figS5_*`, `fig6_*`).

Layout: row A keeps **three species columns** — Arabidopsis (three biological
replicates) | Mouse | HeLa (three replicates) — because it plots curves against
*k*.  The mouse column carries **both independent studies in the same axes and
they are never averaged**: study A = SRP166020 (solid lines, filled dots,
lighter shade), study B = SRP357195 (dashed lines, open dots, darker shade);
every dot and mean bar is drawn once per study, and each study keeps its own
selected combination (a mouse study is a single sequencing unit, so it
contributes one dot per position).

Rows B, C and D are **three full-width tiers sharing one x axis** — the 13 m⁶A
tool configurations, whose names are printed once, **slanted 45°**, right-aligned
at their own tick under row D (no axis title: the names *are* the axis text).
Slanted, neighbouring names are parallel lines ≈ 0.47 in apart — far more than
their cap height — so a name can never collide with its neighbour however long
it is.
Inside every tool position the four independent groups sit as **separate
clusters** (Arabidopsis | mouse study A | mouse study B | HeLa): one dot per
sequencing unit plus the mean *within* that group, and never a mean across
groups or across the two mouse studies.  Row D carries no group dimension at
all: the unmodified controls are not species-specific, so repeating the same
per-tool medians in columns would suggest measurements that do not exist.  Row E
is likewise one full-width panel — the five site sets (2.2 in apart) with their
names in two horizontal lines — the tool names of row D are the only slanted
text on the page.

**Complementary to main-text Figure 6 by construction.**  The rebuilt main
figure carries the selection criterion (Fig. 6D), the selected trade-off over
*k* = 1–5 (Fig. 6A), its per-unit trajectories (Fig. 6B), the configurations'
own recall/PPV and their main effect at fixed *k* (Fig. 6C), the coverage CDF of
the site sets (Fig. 6E) and the control burden of the selected sets (Fig. 6F);
none of those panels or quantities is repeated here.  Figure S5 shows what the
main figure cannot hold: the greedy path beyond five tools, the site quality of
every single configuration, and the DRACH composition of the site sets.

---

**Figure S5. Beyond the five selected tools: the full greedy path, the
per-configuration site quality and the DRACH composition of the site sets.**
Tool combinations were scored on the shared measurable universe of each group
(annotated exons, reference-base compatible, coverage ≥ 10× in every independent
sequencing unit) with a 2-bp GLORI matching window; at each *k* = 1…5 the
combination of the 13 m⁶A tool configurations with the highest mean recall was
selected subject to a mean PPV above the chance level *p*₀ = |GLORI ∩ universe| /
|universe| (1.5 % Arabidopsis; 0.7 % for each mouse study; 0.85 % HeLa).  The
criterion itself — every enumerated combination in the mean-recall versus
mean-PPV plane, *p*₀ and the selected optimum of each *k* — is the subject of
main-text Fig. 6D, and its numbers are the frozen table
`figS5_search_space.tsv`: **all 2,379 enumerated combinations of every group, at
every *k*, satisfy PPV ≥ *p*₀**, so *p*₀ acts as a formal guardrail and the
criterion reduces to the highest mean recall at each *k*; the greedy forward
selection drawn in panel A reproduces all five exhaustive optima (compared as
member sets).  The site-level tables below are read from the exports of
`callsets`, the same clean call layer used by the other supplementary
figures.

**(A)** The full greedy path from the best single tool to all 13
configurations — the only place where *k* > 5 appears (the shaded band marks the
*k* range of the main figure).  Mean union recall of the measurable universe
(filled circles), mean union PPV (open circles), mean intersection PPV (open
squares) and the recall that intersecting costs (open triangles); the dotted
grey line is the best single tool.  Beyond the selected five tools the recall
gain flattens — the remaining eight configurations add +5.7 pp (Arabidopsis),
+4.5 pp (mouse study A), +7.7 pp (study B) and +5.4 pp (HeLa) — while the union
PPV keeps falling (14.1 → 12.6, 6.8 → 6.0, 6.3 → 5.8 and 7.5 → 5.5 %).  The
intersection retains 0.12–2.24 % of the testable GLORI sites at *k* = 5 and
exactly zero at *k* = 13 in every group; its PPV is undefined once the
intersection is empty (no open square after *k* = 7 in Arabidopsis, *k* = 11 in
HeLa and *k* = 12 in study A), which is the quantitative reason why the
intersection strategy is used at the small *k* of the main figure and not
beyond.

**(B)** Coverage (reads, logarithmic) of each configuration's own calls, one dot
per independent unit and a bar for the mean within a group: every
configuration's call set is supported by sequencing depth (group medians 35–863
reads; the least covered configuration still reaches 35 reads in HeLa and 38 in
Arabidopsis).

**(C)** DRACH fraction of the same calls: the configurations separate sharply
(6.0–99.9 % Arabidopsis, 8.2–99.9 % study A, 7.8–99.9 % study B, 7.1–99.9 %
HeLa).  This spread is the site-level basis of the trade-off that panel E
quantifies for the sets the criterion actually selects.

**(D)** Control burden of every configuration on the two unmodified controls —
FP per 10 kb (logarithmic), one dot per sample (filled, Curlcake IVT; open, HeLa
IVT).  Curlcake is the harsher null (median 14.3 versus 0.7 FP/10 kb) and the
three HeLa-only configurations (CHEUI_m6A, DRUMMER, ELIGOS2_diff) have no
Curlcake run, so Curlcake values exist for 10 of the 13 configurations (two of
them as a single sample), whereas HeLa IVT covers all 13 with three replicates.
This is the per-configuration counterpart of main Fig. 6F, which shows the same
controls only for the five selected sets.

**(E)** DRACH composition of the site sets behind the trade-off, one dot per
independent unit and a bar for the mean within a group: the best single tool
(9.5 % Arabidopsis, 99.8–99.9 % elsewhere), the sites the union adds beyond it
(67.3, 37.2, 33.4 and 14.5 %), both intersections (99.97–100 %) and the GLORI
reference (38.8, 90.6, 90.6 and 89.2 %).  The sets differ by three orders of
magnitude in size (marginal sets 132,000–265,000 sites per group; five-tool
intersections 142–593 sites), so the panel reads: the union buys recall with
additions that are predominantly non-DRACH in mouse and HeLa and partially
DRACH in Arabidopsis, while the intersection is a motif-pure, deeply covered
core — the site-level reading of the sensitivity–precision trade-off discussed
in the main text.

---

## Values quoted in the legend (source tables)

| quantity | Arabidopsis | mouse A | mouse B | HeLa | table |
|---|---|---|---|---|---|
| union recall, *k* = 5 → max(*k* = 6–13) | 46.7 → 52.4 % | 80.1 → 84.6 % | 74.0 → 81.7 % | 56.8 → 62.2 % | `figS5_greedy_1to13.tsv` |
| union PPV, *k* = 5 → 13 | 14.1 → 12.6 % | 6.8 → 6.0 % | 6.3 → 5.8 % | 7.5 → 5.5 % | same |
| intersection recall, *k* = 5 → 13 | 0.12 → 0 % | 2.24 → 0 % | 1.74 → 0 % | 1.58 → 0 % | same |
| intersection PPV undefined from *k* = | 8 | 13 | never empty | 12 | same |
| per-tool coverage median, min–max | 38–515 | 58–428 | 44–174 | 35–863 | `figS5_tool_quality.tsv` |
| per-tool DRACH, min–max | 6.0–99.9 % | 8.2–99.9 % | 7.8–99.9 % | 7.1–99.9 % | same |
| FP/10 kb on the controls, median | Curlcake 14.3 / HeLa IVT 0.7 (all groups; controls are not species-specific) | | | | `fig6_negative_control_fp_bytool.tsv` |
| site-set DRACH, single / marginal / isect *k* = 2 / isect *k* = 5 / GLORI | 9.5 / 67.3 / 100 / 100 / 38.8 % | 99.9 / 37.2 / 99.97 / 100 / 90.6 % | 99.9 / 33.4 / 100 / 100 / 90.6 % | 99.8 / 14.5 / 100 / 100 / 89.2 % | `figS5_site_quality.tsv` |
| set size, marginal / isect *k* = 5 (sites) | 264,822 / 142 | 131,515 / 438 | 131,771 / 294 | 168,287 / 593 | same |
| chance level *p*₀ | 1.5 % | 0.7 % | 0.7 % | 0.85 % | `figS5_search_space.tsv` |

Installed as `$RNAMODBENCH_LOCAL/manuscript/sup/sup5.pdf`; the published copy under
`$RNAMODBENCH_LOCAL/submission/` is replaced only on request.
