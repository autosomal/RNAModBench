# Figure 7 (revised) — legend draft (A–D)

**Figure 7. Performance of non-m6A modification detection tools,**
rebuilt at per-replicate resolution from the cleaned site layer
(`harmonisation/callsets`). All metrics are computed per independent
sequencing unit (Human WT n = 3 and unmodified IVT n = 3 units, from two
different studies with the three IVT units sharing one BioSample; independent Curlcake
IVT constructs n = 2 plus one depth-matched subset) inside the candidate-
site universe (coverage >= 10, reference-base compatible). No m6A
reference (GLORI) enters any conclusion in this figure (R2-2).

**(A) Calls per replicate in Human WT versus Human IVT.**

anchored enrichment over chance (right).** Left: in-universe calls of each
replicate (filled, WT; open, Human IVT; short bar, condition mean),
grouped by the three classes the analysis resolves (FP-dominated,
intermediate, specific but sparse). CHEUI-m5C reports 15,845 ± 3,827 calls
per unit in Human WT versus 17,129 ± 7,910 on the Human IVT libraries (mean ± SD
of three units); the raw unions of the three replicates, 47,747 and 51,171,
are kept for reference only. Once the coverage composition of the two libraries
is matched the IVT/WT ratio is 0.79 (95% CI 0.49-1.14), so the union excess
does not survive depth matching.

**(B) Enrichment over chance.** Overlap with external references
(circles, RMBase + DirectRMDB compilation; squares, orthogonal NGS gold standards;
R = 1000, +/-1 bp) divided by its chromosome-stratified permutation
expectation; dashed line, 1 = indistinguishable from random candidate
sites. NanoPsu and NanoSPA-Psi are strongly enriched in WT (mean 84x and
90x; empirical p <= 0.001) with zero or one IVT overlap, whereas CHEUI-m5C
(3.3x vs 2.9x) and NanoMUD-Psi (2.6x vs 2.0x) are equally enriched in WT
and Human IVT. NanoMUD-m1Psi has no external reference (RMBase + DirectRMDB covers
pseudoU, m5C and Nm only) and is therefore absent from this sub-panel.

**(C) False positives on unmodified controls (Curlcake).** Calls per 10^6
candidate sites on the unmodified Curlcake constructs (dots, individual
constructs; open dot, the depth-matched subset of rep3, excluded
from the means; bar, mean of independent constructs). NanoMUD-m1Psi
18,100/18,899, NanoNm 3,903/3,598, NanoMUD-Psi 0/610 and NanoPsu /
NanoSPA-Psi 0/203 per 10^6 candidates. These constructs are covered far more
deeply than the HeLa libraries (median candidate-site coverage 44,755 versus
23), so the panel reads each tool's coverage-dependent detection floor rather
than a uniquely stringent control. CHEUI-m5C made no calls on the Curlcake constructs and is therefore not
shown in this panel.

**(D) Third-party corroboration.** CHEUI probabilities of the
authors' own third-party E. coli data set (GSE271571) for the m5C and m6A
models, E. coli WT versus E. coli IVT (solid, m5C; dashed, m6A); the
high-confidence calls persist on unmodified RNA, i.e. they are a property
of the CHEUI model rather than of our pipeline.

**(E) Replicate consistency.** Within-condition replicate overlap: mean pairwise
Jaccard between the three units of each condition, Human WT (x axis) against the
unmodified IVT libraries (y axis), one circle per tool with whiskers spanning
that condition's three within-unit pairs (log-log; dotted line, WT = IVT).
CHEUI-m5C falls an order of
magnitude below every other tool (mean pairwise Jaccard 0.015 in WT and 0.013 in
IVT versus 0.09-0.41).

**(F) Score validity.**
Mann-Whitney AUC of each tool's own reported score between WT and
Human IVT, with whiskers spanning the nine unit pairs (dotted line,
0.5 = no discrimination). CHEUI-m5C is the
only tool that separates the conditions, and does so in the wrong
direction (AUC 0.21); the Psi and Nm tools stay at 0.50-0.56.

**(G) Replicate-aware metagene distribution of the non-m6A calls**
(1 kb - 5'UTR - CDS - 3'UTR - 1 kb; Ensembl GRCh38.112 annotation;
strand-aware; drawn with the Bioconductor Guitar package on the
replicate-level call sets, as in the original figure). Thick line,
majority consensus (sites present in >= 2 of 3 replicates); thin dashed
lines, individual replicates; blue, Human WT; orange, Human IVT; dotted
verticals, segment boundaries; each panel carries its own
`<tool>-WT` / `<tool>-IVT` key below the axis. Majority site
counts (WT/IVT): CHEUI-m5C 909/776, NanoMUD-Psi 2,566/3,145,
NanoMUD-m1Psi 9,020/11,366, NanoNm 874/2,037, NanoPsu 40/45 and
NanoSPA-Psi 39/46; every panel carries its own key. For every tool the WT
and unmodified-IVT profiles are similar, i.e. the positional profile
carries modification-specific information only where the calls do.
