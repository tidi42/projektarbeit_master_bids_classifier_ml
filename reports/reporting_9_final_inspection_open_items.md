# Reporting 9 — Final pre-submission inspection: open items

Senior-reviewer inspection of **26082026_Projektarbeit.docx** and **26082026_Supplement.docx**
(built from the `_v2` files). This report lists (A) the fixes already applied directly to the two
new documents and (B) the open items that still need author input, a decision, or data I cannot
verify or generate. Items are ordered by severity.

Date: 2026-08-26.

---

## Update 2026-08-26 — author decisions implemented (round 3)

Applied directly to `26082026_Projektarbeit.docx` / `26082026_Supplement.docx`:

- **B2 done** — added a main-text **Table 1** (14-row hypothesis registry, ID / Hypothesis / Family)
  after the *Pre-selected hypotheses* Methods paragraph; rendered and verified.
- **B3 done** — manuscript now names the Kraken2 DB **k2_NCBI_reference_20251007** (was "PrackenDB").
- **B4 resolved → 160 is correct.** The code (`run_mixed_supplement`) fits the model on the
  per-(donor × titration) rows of both compared methods, with and without truncated reads;
  `mixed_model_supplement.tsv` records `n_units = 160` (= 8 × 5 × 2 × 2). The Discussion's "~160" was
  right; the Supplement's "~40" was corrected to 160 with its composition.
- **B5 verified → "longest-running", not "deepest".** BLAST wall-clock (file mtimes): D06 c1 = 233 min
  (3.9 h, the longest) vs D08 c1 = 230 min; but D08 c1 is the read-depth peak (771,919 vs 545,504
  non-human reads). Reworded to "the longest-running library, D06 c1".
- **B6 verified → 2.5% is correct.** `model_comparison.tsv` gives fixed-threshold precision at c5 =
  0.0250. Supplement Fig S3 legend changed "~3%" → "~2.5%" (Abstract already 2.5%).
- **B7 verified → 116,323,055 is correct** (already in the manuscript). Pipeline log: "Stage 05
  complete: 116323055 predictions"; the `.gz` has a header row, so a line count of 116,323,056 is the
  header + 116,323,055 data rows. No change needed.
- **B8 done** — LOEO standardised to *leave-one-experiment-out* throughout both documents.
- **B9 done** — single-sex limitation sentence added to the Discussion.
- **B10 done** — standardised on **NC** (figures already use NC, so no re-render); "c0" text updated.
- **C1 done** — Figure 2A y-axis floor is now data-driven (0.85 here) so the GLM/Kraken2 whisker
  (down to ~0.881) is no longer clipped; the fixed-run figure was regenerated and re-embedded in the
  manuscript, and the calculated-run figure regenerated too.

**Still open: B1** — the reference list remains intentionally empty for now (per author).

---

## A. Fixes already applied directly (no action needed)

Manuscript:
- **Figure 1B** was uncited in the text; added a citation in the Dataset-parameters paragraph, in
  panel order (A → B → C).
- **"no-template control (NC)" → "no-spike control (NC)"** in the Results text and the Figure 1
  legend. The barcode is unspiked donor blood (contains donor DNA/microbiome/kitome), so
  "no-template" was a misnomer; this also aligns with the "no-spike control(s)" wording already used
  elsewhere in Results.
- **"the hac basecaller" → "the HAC basecaller"** (Discussion) — casing consistent with the defined
  abbreviation.
- **"per- species" → "per-species"** (Discussion) — stray space.
- Reference-preparation parenthetical **"(17 species total -> zymo_members.fasta)" →
  "(17 species-level taxa total)"** — removed the informal code-style arrow/filename from main-text
  prose.

Supplement:
- **Factual correction (c1/c5):** "pooled metrics would otherwise be dominated by the highest-depth
  **c5** libraries" → **c1**. c1 is the highest-input level and holds 61.8 % of the reads; c5 is the
  lowest-input, read-poorest level.
- **Stale contradiction fixed (GLMM training set):** the Cross-validation-design paragraph said the
  mixed-effects model "was accordingly fit on the complete, undownsampled read set", contradicting the
  manuscript (Methods and Discussion), which correctly state the GLMM used a **250 000-read
  label-stratified subsample**. Rewritten to: trees/GLM/threshold on the complete set; the GLMM alone
  on a label-stratified subsample of ≤ 250 000 reads/fold (natural class balance preserved, no
  prior-probability correction).
- **Broken section cross-reference:** "(Pre-registered hypotheses, H9b)" → "(Pre-selected
  hypotheses, H9b)" — the section is titled *Pre-selected hypotheses*.
- **Dangling item reference:** "(see Caveats and limitations, item 2)" → "(see Caveats and
  limitations)" — the Caveats section is not itemised/numbered.

All edits preserved the citation/field codes and all embedded figures; both documents are valid and
render. **After opening in Word, refresh both Tables of Contents (Ctrl+A → F9)** because content
lengths changed (also listed in B-12).

---

## B. Open items still to do (author input / decision / data required)

### Blockers — must be resolved before submission

**B-1. Reference list is empty.**
The *References* section contains no entries and there are no in-text citations anywhere in the
manuscript. Add citations and a reference list covering at least: minimap2, BLAST+/BLASTn, Kraken2,
ranger, XGBoost, glmmTMB, lmerTest, seqkit, samtools, barrnap, dustmasker, Dorado; the ZymoBIOMICS
D6331 standard; the GRCh38 and T2T-CHM13 references; NCBI core_nt; and the statistical methods
(Wilcoxon signed-rank, Holm–Šídák, Platt scaling, Poisson detectability, ONT adaptive sampling).
This is the single largest gap. (I cannot generate citations without fabricating them.)

### Style & completeness — recommended

**Terminology — "pre-selected" vs "pre-registered".** The two terms are used interchangeably across
the documents; standardise on "pre-selected" unless a formal pre-registration was filed.

**B-11. Figure/table image content.** Figure 2A's clipped y-axis is now fixed and re-embedded; before
submission confirm the remaining figures match their legends — correct panel labels, **species names
in italics** (see B-13), no stray artefacts, and axis units matching the text.

### Housekeeping

**B-12. Refresh the Tables of Contents.** Open each document in Word and update the TOC
(Ctrl+A → F9); page numbers shifted with the edits.

**B-13. Italicise species names** (optional but expected by most journals): *Escherichia coli*,
*Candida albicans*, *Saccharomyces cerevisiae*, etc., throughout both documents and in figures.

---

*Numbers verified as internally consistent during this pass (no action needed):* the read-accounting
chain (66,058,488 sequenced − 61,310,642 human = 4,747,846 non-human; − 24,864 ambiguous/indeterminate
= 4,722,982 combined-arm; 57 % positive), c1 = 2,932,304 (61.8 %), the 52-fold donor×level spread, all
per-hypothesis effect sizes and p-values across Results / Table 2 / Supplementary Tables S2–S3, the
baseline gain of 0.054, and the Poisson E[N] ≥ 3 at P ≥ 0.95.
