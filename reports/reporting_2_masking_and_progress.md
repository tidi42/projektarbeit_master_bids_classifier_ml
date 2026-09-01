# Reporting 2 — Ambiguous-region masking upgrade + pipeline progress reporting

**Addendum to** `reports/reporting_1_methodology.md`. This report documents the
changes made after reporting 1:

- **(A)** the ground-truth exclusion BED was **rebuilt and extended** — most
  importantly with **dustmasker low-complexity masking** (adopted at **high
  priority**; reasoning in §3), plus **plasmid-backbone masking**;
- **(B)** a **regression fix** to the BED introduced by the genome collapse;
- **(C)** **pipeline progress reporting** so a long run shows which stage is
  running and which are done.

It supersedes reporting 1 §5 (last paragraph), §19 items 8 and 9 (see §7).
Document date: 2026-08-03. Companion trackers: `open_items_v3.md`.

---

## 1. Summary of changes

| Change | File(s) | Effect |
|---|---|---|
| Rebuild BED against the **collapsed, species-named** genomes | `prepare_ambiguous_regions.R` | BED coordinates now match `ref/genomes/<species>.fasta`; fixes the regression (§2). |
| **Dustmasker low-complexity** masking, ≥ 300 bp blocks | `prepare_ambiguous_regions.R` → `dustmasker_lowcomplexity()` | Substantial low-complexity / simple-repeat loci excluded from the positive ground truth (§3). |
| **Plasmid-backbone** masking | `prepare_ambiguous_regions.R` → `plasmid_regions()` | Plasmid contigs (lost to the genome collapse) re-masked in collapsed coordinates (§4). |
| barrnap rRNA operons (already in reporting 1) | `prepare_ambiguous_regions.R` → `barrnap_operons()` | Full rRNA operons unioned into the BED. |
| New parameter `lowcomplexity_min_len = 300L` | `00_config.R` | Region-mask threshold that separates substantial repeats from short tracts (§3). |
| Open item 2 → **resolved** | `00_config.R` (`cfg$open_items`) | Only mobile-element masking now deferred (optional). |
| **Progress reporting** (plan, per-stage bar + timing, completion summary; per-library bars) | `run_pipeline.R`, `utils.R`, `02/03/05_*.R` | The run shows which process is working and which are completed (§5). |

**Resulting BED:** `ref/zymo_ambiguous_regions.bed` = **124 merged, species-named
regions** = union of 16S/18S gene bodies + barrnap operons + plasmid backbones +
dustmasker low-complexity (≥ 300 bp).

---

## 2. Why the BED had to be rebuilt (regression fix)

`prepare_zymo_genomes.R` rebuilds the reference as **collapsed, canonical-species**
FASTAs (the 5 E. coli strains → one `Escherichia_coli`; supplier headers normalised).
The previous BED had been written against the **raw** contig names
(`E.coli.JM109.final.genome`, the misspelled `Clostridium_difficille`, …), so after the
collapse those names no longer existed in `ref/genomes/*` and the BED could not
intersect the reads — a silent regression that would have disabled all region masking.

**Fix:** `prepare_ambiguous_regions.R` now derives regions **against the built
`ref/genomes/<species>.fasta`** (sourcing `canonical_species` from
`prepare_zymo_genomes.R`), so BED names are exactly the species contigs used at
alignment time.

**Verification:** the 124 region names were diffed against the FASTA members with
`comm -23` → **empty** (every region maps to a real member; **0 stale names**).

---

## 3. Dustmasker low-complexity masking — adopted at HIGH priority

**What it does.** `dustmasker` (BLAST+ SDUST) flags low-complexity / simple-repeat
tracts — homopolymer runs, di-/tri-nucleotide microsatellites, satellite-like arrays —
in each Zymo genome. Substantial blocks (**≥ 300 bp**, `cfg$params$lowcomplexity_min_len`)
are unioned into the exclusion BED, so any read whose **best Zymo hit falls inside such a
tract is excluded from the ground-truth positive set** (it is not forced into the binary).

**Why this is high priority — the core argument.**

1. **Low-complexity ≈ non-specific.** Simple-repeat tracts are shared across genomes and
   are compositionally generic. A human or donor-bacterial read sitting over a
   low-complexity stretch can align to a Zymo genome at **high nominal identity by
   composition, not by true homology**. In the ground truth that alignment would be
   scored as a **spurious Zymo positive**.

2. **It would make the ground truth circular.** The entire project is about separating
   **true Zymo reads from false positives** in a ~99 % human, low-biomass matrix. If the
   positive **label** itself admits low-complexity artefact alignments, we would be baking
   the very false positives we are trying to detect **into the target** — inflating the
   apparent positive set and corrupting every downstream number (precision@recall, FDR)
   and every hypothesis test that compares classifiers against that target.

3. **It hits exactly the metric we care about most.** The clinical operating points are
   **high recall (0.95 / 0.99)**; precision and FDR there are dominated by the *tail of
   borderline positives* — which is precisely where low-complexity cross-matches live.
   Leaving them in would most damage the numbers the Projektarbeit is judged on.

4. **It ties directly to the contamination hypotheses.** Low-complexity cross-matching is
   a textbook source of taxonomic mis-assignment (relevant to H7 and the core_nt
   contamination angle). Removing these loci from the truth keeps the read-level positive
   definition defensible and keeps the false-positive analysis about *real* contamination
   rather than alignment artefacts.

⇒ In a 99 %-human, low-biomass regime evaluated at high recall, unmasked low-complexity
artefacts attack the **integrity of the positive class** in exactly the regime that
defines the project. That is why low-complexity masking is treated as **high priority**,
not an optional refinement.

**Why a ≥ 300 bp region-mask threshold (not "mask everything").** Raw dustmasker output on
the Zymo genomes is **~99,774 micro-intervals with a mean length ≈ 37 bp**. Excluding a
read on *any* overlap with *any* of these would clip **nearly every long ONT read** (a
multi-kb read spans many short tracts) and would discard almost the whole positive set —
destroying recall. So the design is **two-tier**:

- **Substantial loci (≥ 300 bp)** — satellites, long homopolymer / SSR arrays — are
  **region-masked** (excluded from ground truth). This keeps **107 substantial blocks**,
  merged into the final BED.
- **Short, ubiquitous stretches (< 300 bp)** are **not** region-masked; they are handled
  **per read** by the complexity **features** the model already carries
  (`dust_score`, `shannon_entropy`, `homopolymer_frac`, GC). The model can learn to
  down-weight a read with a short low-complexity segment **without throwing the read away**.

This split removes the circular, artefact-inducing loci **while preserving recall** — the
big repeats are excluded from truth, the small ones are left to the classifier.

---

## 4. Plasmid-backbone masking (recovered after the collapse)

The genome collapse merges each species' contigs (chromosome + plasmids) into one FASTA
with N-spacers, which **loses the plasmid headers** — so header-based plasmid masking was
no longer possible (this was the "plasmid deferred" caveat in reporting 1). `plasmid_regions()`
now recomputes each plasmid's **offset in the collapsed coordinate system** (accounting for
the 100 bp N-spacer between concatenated contigs) and masks those spans.

**Reasoning:** plasmids are frequently **horizontally shared and high-copy**, so they drive
cross-species alignments and inflate positives the same way rRNA does; keeping them out of
the read-level positive truth is consistent with the rest of the ambiguity policy.

---

## 5. Pipeline progress reporting (UX)

A long run (external tools over 48 libraries, then training) previously gave little live
feedback. `run_pipeline.R` now reports **which stage is running and which are done**:

- **Run plan** up front — the ordered list of stages that will execute.
- **Per-stage banner** — a progress bar plus `RUNNING` on entry and `DONE [duration]` on
  exit, so overall pipeline position is always visible.
- **Completion breakdown** — a `[x]` line per stage with its wall-clock time and a total.
- **Within stages** — `02_ground_truth_labels.R` and `03_build_features.R` show a
  **per-library** `txtProgressBar` (via the new `make_progress()` in `utils.R`);
  `05_train_models.R` prints a per-arm training message.

Illustrative console output:

```
-- RUN PLAN: 7 stage(s) --------------------------------------------
   1/7  stage 01  external tools (QC/minimap2/BLAST/Kraken2)
   ...
   7/7  stage 07  hypothesis tests H1-H12
--------------------------------------------------------------------

[#######...............]  14%  (1/7)  RUNNING  stage 01 : external tools ...
[#######...............]  14%  (1/7)  DONE     stage 01 : external tools ...  [42.0m]
...
==================== PIPELINE COMPLETE ====================
   [x] stage 01  external tools ...                      42.0m
   ...
   total: 1.35h
==========================================================
```

This is purely observability — no analysis logic changed. `make_progress(total, label)`
is a no-op when `total <= 0` and degrades cleanly under `Rscript`.

---

## 6. Verification

- **Syntax:** all edited scripts pass `Rscript -e "invisible(parse(...))"`
  (`00_config.R`, `utils.R`, `run_pipeline.R`, `02/03/05_*.R`).
- **BED consistency:** 124 merged regions; every region name ∈ FASTA members
  (`comm -23` empty); 0 stale contig names; coordinates within contig lengths.
- **Progress visuals:** stage bar, run plan and `make_progress()` bar smoke-tested.

---

## 7. Supersedes in reporting 1

- **§5 (last paragraph)** — "*Deferred: plasmid backbones and mobile elements*" →
  plasmid backbones are **now masked**; only **mobile elements** remain deferred (optional).
- **§19 item 8** — plasmid masking is **done**; only mobile-element masking is deferred.
- **§19 item 9** — `ref/zymo_members.fasta` is **built** and the **Kraken2 DB is now
  provided** (`k2_NCBI_reference_20251007`, extracted under `database_kraken2/`, wired to
  `cfg$paths$kraken2_db`); **no hard full-run prerequisite remains** (`open_items_v3.md` §1).
- **`cfg$open_items` item 2** — now **resolved**.

---

## 8. Reproduce

```bash
Rscript scripts/prepare_zymo_genomes.R          # collapsed, species-named genomes + contig2species
Rscript scripts/prepare_ambiguous_regions.R     # 16S/18S + barrnap + plasmid + dustmasker(>=300bp) BED
# delete ref/zymo_ambiguous_regions.bed to force a full rebuild
```

Parameters live in `00_config.R` (`lowcomplexity_min_len = 300L`); the reference genomes
must be built first because the BED is derived against `ref/genomes/*`.

---

## 9. Scope notes (2026-08-03)

**Note — H7 (adapter content) deferred to a later experiment.** H7 asked whether reads
carrying residual adapter/barcode sequence are enriched among false positives (Fisher / odds
ratio + FDR ablation). It is a **separate experiment for later** and is **not run in the
current pipeline**: `cfg$hypotheses` marks H7 `active = FALSE`, `07_hypothesis_tests.R` emits
a *DEFERRED* provenance row only (mirroring the dropped H8), and the per-read
`work/adapter_content.tsv` input is no longer required for now. The primary family (H1–H6,
alpha-spent) and the remaining secondary tests (H9–H12) are unaffected. *Re-enable* by
restoring the Fisher/OR body in `test_H7`, setting H7 `active = TRUE`, and producing the
`adapter_positive` table. Rationale: adapter content is orthogonal to the core
Zymo-vs-false-positive discrimination and belongs with the core_nt contamination follow-up,
which is scoped as future work.

**Full-run analysis only (no subsampling).** The pipeline runs end-to-end on the **complete**
human-depleted read set; there is **no development subsample**. Open items 11/12 (subsample
feasibility / acceptability) are resolved as **full-run only** — BLASTn, Kraken2 and all
downstream metrics score every non-human read. This is consistent with the no-downsampling
stance for the negative class (`negative_downsample_ratio = NA`): at pathogen-relevant
abundance the signal is a few reads per million, so discarding reads would remove the very
events the classifier must catch.
