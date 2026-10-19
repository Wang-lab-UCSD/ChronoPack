# Time-tag-bias filter: bulk vs single-cell negative control

Status (2026-10-10): **observations only, no change made.** The pipeline still flags
time-tag-biased genes against the bulk control (`analyses/inputs/k562_neg-ctrl_mtx-agg.pkl`).
Changing the control (or combining controls) means re-running from `03_PrepareForTradeSeq`
onward.

Analysis files: `verification_runs/negctrl_provenance_20261010/` (bulk control provenance) and
`verification_runs/negctrl_sc_vs_bulk_20261010/` (comparison; per-gene tables in `outputs/`).

## The two controls

| | Bulk control (current) | Single-cell control |
|---|---|---|
| Source | N706 "Bulk K562 Control 600 Cells", `~/projects/chrono_seq/analyses/old/k562-TNFa.hg38/N706_Bulk_K562_Control_600_Cells.aligned.sorted.merged.tagged.dge.time_tags.csv.gz` | `~/projects/chrono_seq/analyses/qc/N701-N704_scK562_NegCTRL_Rev_Time_2X_Cells_RT{1-4}...time_tags.csv.gz` (tags added in reverse time order) |
| Pseudobulk | all 35,030 barcodes summed per `FINAL TIME TAG` | all 9,387 barcodes of the 4 libraries summed per `FINAL TIME TAG` |
| Counts per tag | 1.0-1.5M | 2.0-3.8M |
| Genes | 27,085 | 30,501 (23,929 in both) |

The bulk pickle was written by commented-out code in chrono_seq
`thp1-k562-lps.hg38.rep/notebooks/PrepareForTradeSeq.ipynb` (cells 4-5) and rebuilds
byte-identically from the N706 file. The single-cell pseudobulk follows chrono_seq
`k562-tnfa.hg38.rep/notebooks/not_used/old/BetterInterventionDesign.ipynb` (`rev_neg_ctrl`).

## Observations

- **Pairing.** `FINAL TIME TAG` behaves as the tag identity: the two controls' per-gene tag
  profiles agree better at the same tag index (median r = 0.14) than reversed (0.02), and the
  single-cell control flags more genes with the same index than reversed (K562a 1,379 vs 523).
- **Agreement.** The controls agree only partly. On 19,007 K562a genes scored by both, 313 are
  flagged by both versus ~72 expected by chance (4.3x). Among genes that would reach tradeSeq:

  | Context | Excluded, bulk | Excluded, single-cell | Both | Jaccard | Single-cell would drop (r9-significant) | Would return |
  |---|---|---|---|---|---|---|
  | K562a | 753 | 1,075 | 253 | 0.16 | 822 (219) | 500 |
  | K562b | 563 | 838 | 173 | 0.14 | 665 (127) | 390 |
  | THP1 | 262 | 309 | 38 | 0.07 | 271 (13) | 224 |

- **Shared core.** Genes flagged by both are mostly abundant housekeeping transcripts (hnRNPs,
  HSPA5/8/9, CANX, PDIA6, translation factors, proteasome subunits, YBX1, PTMA, FTL, HMGB1) plus
  some processed/mitochondrial pseudogenes; list: `outputs/excluded_by_both.K562a.csv`.
- **Biologically expected genes removed by one control only (K562a):**
  - bulk only: erythroid program HBG1, HBG2, GYPA (HBG2 is the erythroid-like fate marker);
  - single-cell only: TFRC; MXD1 (meanLogFC 0.70, Wald 62).
  - Flagged by both: STAT1, RIPK1, MAPKAPK2, TGFB1, VEGFA; SPTA1, SPTB, ANK1, FTL.
  - No classic fast NF-kB target (NFKBIA, TNFAIP3, ICAM1, CXCL8, BIRC3, IRF1, ...) is flagged
    by either control.
- **Strength of control-specific exclusions (K562a).** Single-cell-only exclusions include 219
  r9-significant genes, several very strong (EIF3B Wald 219.5, HK1, TFRC, NUDT4); bulk-only
  exclusions include 60 genes with p < 0.05 (older K562a table), mostly weaker (H2AC12 the
  exception). Ranked lists: `outputs/{bulk_only,sc_only}.K562a.ranked.csv`.
- **False-positive estimate.** The reversed pairing of the single-cell control is a mismatched
  baseline: it flags ~40% as many genes as the correct pairing, so a similar share of
  single-cell flags may be chance at the current thresholds (r >= 0.7, explained >= 0.25). No
  such estimate exists for the bulk control yet.

## Options considered (not decided)

1. Keep the bulk control (current).
2. Exclude genes flagged by **both** controls (K562a 253 / K562b 173 / THP1 38): keeps the
   shared housekeeping core out, retains HBG1/HBG2/GYPA, MXD1, TFRC. Leading candidate.
3. Single-cell control only: better matched (same single-cell workflow and cell line, deeper,
   4 libraries, built-in false-positive estimate) but more aggressive.
4. Union of both: strictest.

Possible next checks: mismatched-pairing false-positive rate for the bulk control; tag profiles
of HBG2, MXD1 and TFRC under both controls.

## Related fix

`load_negative_control_logcpm` now casts the control table to one float block and
`compute_smooth_correlations` aligns the tables once before the gene loop. The single-cell
pseudobulk (summed from libraries with different genes) had 29,825 mixed int/float blocks,
which made one pass take ~3.6 h instead of ~2 min. Results unchanged.

## Shared time trends among K562a significant genes (2026-10-10)

The filter above uses tag-to-tag (residual) correlations only. For the 2,022 K562a genes with
screen p < 0.05, the spline-smoothed tag profiles of the stimulated pseudobulk were correlated
with each control's (`verification_runs/negctrl_sc_vs_bulk_20261010/shared_trends.py`;
tables `outputs/shared_trends.K562a*.csv`).

- Bulk control: 371 genes with trend r >= 0.8 versus 37 with r <= -0.8 (median r 0.42), a clear
  excess of shared positive trends. 299 of them also change no more than 1.5x as much in the
  stimulated cells as in the control (111 with |peak log2FC| >= 0.5), e.g. EGR1, SQSTM1, PHACTR1,
  NFE2L2, IER2, CORO1C, PARK7, ID2, SLC38A2, HK1, TFRC, SOCS1/2, NUDT4.
- Single-cell control: symmetric (204 at r >= 0.8, 234 at r <= -0.8, median -0.01), i.e. about
  what smoothed 12-point profiles give by chance; its trend correlations are uncorrelated with the
  bulk ones across genes (r = 0.0). So it does not confirm or refute the bulk-shared trends.
- The bulk control's tags were added in the same order as the experiments' and the single-cell
  control's in reverse, so a shared bulk trend may reflect handling / time-of-addition effects
  (e.g. immediate-early genes EGR1, IER2, ID1-3) as well as tag effects; this is a hypothesis,
  not tested.
- Canonical NF-kB targets often share the direction of the bulk control's drift (NFKBIA r 0.76,
  CD83 0.90, IER3 0.88, IRF1 0.76) but change 2.8-7x more in the stimulated cells.

## Prototype: stimulation x time interaction vs the negative controls (2026-10-10, exploratory)

`verification_runs/stim_vs_control_interaction_20261010/` (`prepare_pseudobulks.py`,
`interaction_model.R`; tables in `outputs/`, diagnostics and example trajectories in `plots/`).
Library x tag pseudobulks, edgeR QL with `~ 0 + library + ns(elapsed, 4) + stim:ns(elapsed, 4)`;
Delta(t) = [f_stim(t) - f_stim(0)] - [f_ctrl(t) - f_ctrl(0)] from the interaction coefficients, with
simultaneous 90% bands. Elapsed time: stim and bulk control tag t = (t-1)*10 min; single-cell
control (reverse-time) ASSUMED tag t = (12-t)*10 min.

- Replication: each condition is one run (4 RT libraries = technical replicates; bulk control 1
  library), so condition is confounded with run. Within-run fake split: 0 interactions at FDR 0.05,
  2.7% at p < 0.05 (calibrated for technical noise). Two unstimulated runs (single-cell vs bulk
  control): 2.1% at FDR 0.05, 14.7% at p < 0.05, so interaction p-values against a separately run
  control are anti-conservative. Elapsed-time alignment gives fewer spurious interactions than
  tag alignment (155 vs 259), consistent with the assumed reverse mapping.
- Point estimates of max|Delta| are noisy (> 0.5 for 42% of genes in the within-run null); a band rule
  is used instead: response (band outside +-0.5 at some t), control-like (band inside +-0.5 for all
  t), ambiguous.
- Single-cell control, 2,020 nominal genes: 185 response, 227 control-like, 1,608 ambiguous (FDR set
  1,043: 150 / 129 / 764). In the within-run null only ~10% of genes can be shown control-like, so
  "ambiguous" mostly reflects precision. Bulk control: 63 / 89 / 1,864 (single library, wide bands);
  no gene is response under one control and control-like under the other.
- NFKBIA, CD83, IER3, IRF1: response vs the single-cell control (max|Delta| 3.5-4.6 log2, lower bands
  2.4-3.9); vs bulk, CD83 is ambiguous (wide band). The 299 correlation-flagged genes: 45 response,
  38 control-like, 216 ambiguous.

### Equivalence test, +-0.5 log2 with BH (2026-10-12, exploratory)

`equivalence_test.R` (tables `outputs/equivalence.*.csv`, `outputs/equivalence_summary.csv`).
Intersection-union TOST: at every elapsed minute a TOST p for H0 |Delta(t)| >= 0.5; gene p = max
over t; BH within the gene set. A "variance-inflated" sensitivity run scales the covariance by the
between-run lambda (1.29, from ctrl_sc vs ctrl_bulk interaction F; within-run split 0.85).

- Excluded at q < 0.05 (nominal / FDR set / of the 299): single-cell control 168 / 104 / 23
  (inflated 115 / 74 / 16); bulk control 0 / 0 / 0 (smallest p 7e-5 with one library).
- Every gene excluded has a stimulated max change < 0.43 log2 and a control change < 0.5 (129 < 0.25):
  mostly abundant housekeeping transcripts (MALAT1, RPLP0/2, UBB, TUBA1B, FTH1, ATP5*, COX*). Of
  the ~990 genes moving >= 0.5 log2 in the stimulated cells, none reach q < 0.05 (4 have p < 0.05).
  The test therefore removes low-amplitude genes; it does not show that the control explains a trend.
- NFKBIA, CD83, IER3, IRF1: p_equiv ~1 with either control. NFE2L2 q 0.088, PARK7 q 0.067 (sc).

### Normalized trajectory discrepancy D (2026-10-12, exploratory)

`similarity.R` (+ `prepare_bulk_halves.py`); tables `outputs/similarity.*.csv`,
`outputs/similarity_summary.csv`; plots `plots/similarity.D_hist.png`, `plots/similarity.examples.png`.
D = rms_t[s(t) - c(t)] / rms_t[s(t)] (t = 1..110 min, log2 change from t0, same edgeR model). D_noise =
expected D if the true trajectories were identical; D_upper90 from 1,000 coefficient draws (conservative).

- Noise dominates D for most genes: among genes with stimulated change >= 0.5, median D_noise is 1.03
  (bulk control) and 0.73 (single-cell control). Point-D cutoffs therefore mostly select chance low
  values (bulk D < 0.4: 83 genes, e.g. IL6 D 0.08 with D_noise 1.6); D_upper90 < 0.4: 0 (bulk), 1 (sc).
- Precisely measured genes (D_noise < 0.3): bulk 20, sc 48. Low D among them: bulk NFE2L2 0.24, EGR1 0.27;
  sc MIR223HG 0.18, NFE2L2 0.23, EGR1 0.26, ID3 0.27, SOCS1 0.35, NAB2 0.37 (TFRC, ARL4A, PHACTR1, FBXO5,
  ID1 < 0.5).
- Tag-aligned single-cell control reverses the shared trends (EGR1, ID2, ID3, IER2 shape r -0.92 to -0.99;
  D 1.7-2.6): the shared trends follow elapsed time since handling, not tag identity, i.e. they are not
  time-tag bias.
- Named genes (D bulk / sc): EGR1 0.27 / 0.26; ID1 0.30 (noise 0.37) / 0.48; ID2 0.56 / 0.80; ID3 0.71 / 0.27;
  IER2 0.98 / 1.15; NFKBIA 0.89 / 0.94; CD83 0.58 / 1.02 (bulk control rises 1.7 log2); IER3 0.87 / 0.70;
  IRF1 0.94 / 0.98.
- Bulk noise check: half-split (sampling-only) curve noise is 0.62-0.67x the model prediction, which
  borrows dispersion from the stimulated libraries; library-prep noise is not in the half split.
