from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/new-stg/home/alex/tmp/ChronoPack")
RUN = ROOT / "verification_runs" / "gene_dynamics_migration_20260918"
REFERENCE = (
    ROOT
    / "verification_runs"
    / "filename_migration_20260917_074750"
    / "k562-tnfa-puromycin.hg38"
    / "intmd"
    / "sc-K562-1e6-to-tradeSeq.fixed"
)
PACKAGED = RUN / "inputs_packaged_tradeseq_ce"
UNIFIED = RUN / "inputs_unified_chronopack"


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def numeric_summary(left: pd.DataFrame, right: pd.DataFrame) -> dict[str, object]:
    common_rows = left.index.intersection(right.index)
    common_cols = left.columns.intersection(right.columns)
    diffs = {}
    for col in common_cols:
        l = pd.to_numeric(left.loc[common_rows, col], errors="coerce")
        r = pd.to_numeric(right.loc[common_rows, col], errors="coerce")
        both = l.notna() & r.notna()
        if both.any():
            delta = (l[both] - r[both]).abs()
            corr = np.corrcoef(l[both], r[both])[0, 1] if both.sum() > 1 else np.nan
            diffs[col] = {
                "max_abs_diff": float(delta.max()),
                "median_abs_diff": float(delta.median()),
                "correlation": float(corr),
            }
    return {
        "shape_left": list(left.shape),
        "shape_right": list(right.shape),
        "same_index": left.index.tolist() == right.index.tolist(),
        "same_columns": left.columns.tolist() == right.columns.tolist(),
        "row_overlap": int(len(common_rows)),
        "column_overlap": int(len(common_cols)),
        "numeric_columns": diffs,
    }


def rds_dims(path: Path) -> dict[str, object]:
    code = f"""
    suppressPackageStartupMessages({{library(SingleCellExperiment)}})
    x <- readRDS({str(path)!r})
    if (inherits(x, "SingleCellExperiment")) {{
      cat(paste(c(dim(x), length(rownames(x))), collapse="\\t"))
    }} else if (is.list(x)) {{
      cat(paste(c(length(x), length(names(x))), collapse="\\t"))
    }} else {{
      cat(class(x)[1])
    }}
    """
    proc = subprocess.run(
        ["/new-stg/home/alex/miniconda3/envs/tradeseq_ce/bin/Rscript", "-e", code],
        check=True,
        text=True,
        capture_output=True,
    )
    return {"path": str(path), "summary": proc.stdout.strip()}


def rds_name_compare(left: Path, right: Path) -> dict[str, object]:
    code = f"""
    suppressPackageStartupMessages({{library(SingleCellExperiment)}})
    get_names <- function(path) {{
      x <- readRDS(path)
      if (inherits(x, "SingleCellExperiment")) {{
        rownames(x)
      }} else if (is.list(x)) {{
        names(x)
      }} else {{
        character()
      }}
    }}
    a <- get_names({str(left)!r})
    b <- get_names({str(right)!r})
    common <- intersect(a, b)
    cat(paste(
      length(a),
      length(b),
      length(common),
      identical(a, b),
      paste(head(a, 5), collapse = ","),
      sep = "\\t"
    ))
    """
    proc = subprocess.run(
        ["/new-stg/home/alex/miniconda3/envs/tradeseq_ce/bin/Rscript", "-e", code],
        check=True,
        text=True,
        capture_output=True,
    )
    n_left, n_right, n_common, same_order, head_names = proc.stdout.strip().split("\t")
    return {
        "left": str(left),
        "right": str(right),
        "n_left": int(n_left),
        "n_right": int(n_right),
        "n_common": int(n_common),
        "same_order": same_order == "TRUE",
        "head": head_names.split(",") if head_names else [],
    }


def predict_smooth_compare(left: Path, right: Path, *, n_points: int = 9) -> dict[str, object]:
    code = f"""
    suppressPackageStartupMessages({{library(tradeSeq)}})
    as_matrix <- function(x, genes) {{
      if (is.data.frame(x)) {{
        key_cols <- intersect(c("gene", "lineage", "time"), names(x))
        value_cols <- setdiff(names(x), key_cols)
        numeric_cols <- value_cols[vapply(x[value_cols], is.numeric, logical(1))]
        if (!length(numeric_cols)) stop("No numeric columns in tidy predictSmooth output.")
        x$.rowkey <- do.call(paste, c(x[key_cols], sep = "\\r"))
        x <- x[order(x$.rowkey), , drop = FALSE]
        mat <- as.matrix(x[numeric_cols])
        rownames(mat) <- x$.rowkey
        return(mat)
      }}
      if (is.matrix(x)) {{
        return(x[intersect(genes, rownames(x)), , drop = FALSE])
      }}
      if (is.list(x)) {{
        mats <- lapply(seq_along(x), function(i) {{
          m <- x[[i]]
          if (is.null(dim(m))) stop("Unsupported predictSmooth list element.")
          m <- as.matrix(m)
          rownames(m) <- paste0(rownames(m), "\\r", names(x)[i])
          m
        }})
        return(do.call(rbind, mats))
      }}
      stop("Unsupported predictSmooth output class: ", paste(class(x), collapse = ","))
    }}
    g1 <- readRDS({str(left)!r})
    g2 <- readRDS({str(right)!r})
    genes <- intersect(names(g1), names(g2))
    if (!length(genes)) stop("No shared genes in gamList objects.")
    p1 <- as_matrix(tradeSeq::predictSmooth(g1, gene = genes, nPoints = {int(n_points)}), genes)
    p2 <- as_matrix(tradeSeq::predictSmooth(g2, gene = genes, nPoints = {int(n_points)}), genes)
    common_rows <- intersect(rownames(p1), rownames(p2))
    common_cols <- intersect(colnames(p1), colnames(p2))
    if (!length(common_cols)) common_cols <- seq_len(min(ncol(p1), ncol(p2)))
    p1 <- p1[common_rows, common_cols, drop = FALSE]
    p2 <- p2[common_rows, common_cols, drop = FALSE]
    deltas <- abs(as.numeric(p1) - as.numeric(p2))
    vals1 <- as.numeric(p1)
    vals2 <- as.numeric(p2)
    keep <- is.finite(vals1) & is.finite(vals2)
    cat(paste(
      length(genes),
      nrow(p1),
      ncol(p1),
      max(deltas, na.rm = TRUE),
      median(deltas, na.rm = TRUE),
      if (sum(keep) > 1) cor(vals1[keep], vals2[keep]) else NA_real_,
      sep = "\\t"
    ))
    """
    proc = subprocess.run(
        ["/new-stg/home/alex/miniconda3/envs/tradeseq_ce/bin/Rscript", "-e", code],
        check=True,
        text=True,
        capture_output=True,
    )
    n_genes, n_rows, numeric_cols, max_abs, median_abs, corr = proc.stdout.strip().split("\t")
    return {
        "left": str(left),
        "right": str(right),
        "n_genes": int(n_genes),
        "n_rows": int(n_rows),
        "n_columns": int(numeric_cols),
        "max_abs_diff": float(max_abs),
        "median_abs_diff": float(median_abs),
        "correlation": float(corr),
    }


def main() -> None:
    summary = {}

    for filename in [
        "counts_genes_x_cells.csv.gz",
        "cell_times.csv.gz",
        "subset_genes.csv",
        "cell_libsize_full.csv",
    ]:
        ref = REFERENCE / filename
        packaged = PACKAGED / filename
        unified = UNIFIED / filename
        if filename.endswith(".csv.gz") or filename.endswith(".csv"):
            if filename == "counts_genes_x_cells.csv.gz":
                ref_df = pd.read_csv(ref, index_col=0)
                pkg_df = pd.read_csv(packaged, index_col=0)
                uni_df = pd.read_csv(unified, index_col=0)
                summary[f"input_{filename}_packaged"] = {
                    "same_shape": ref_df.shape == pkg_df.shape,
                    "same_index": ref_df.index.tolist() == pkg_df.index.tolist(),
                    "same_columns": ref_df.columns.tolist() == pkg_df.columns.tolist(),
                    "max_abs_diff": float((ref_df - pkg_df).abs().to_numpy().max()),
                }
                summary[f"input_{filename}_unified"] = {
                    "same_shape": ref_df.shape == uni_df.shape,
                    "same_index": ref_df.index.tolist() == uni_df.index.tolist(),
                    "same_columns": ref_df.columns.tolist() == uni_df.columns.tolist(),
                    "max_abs_diff": float((ref_df - uni_df).abs().to_numpy().max()),
                }
            else:
                summary[f"input_{filename}_packaged_identical_bytes"] = ref.read_bytes() == packaged.read_bytes()
                summary[f"input_{filename}_unified_identical_bytes"] = ref.read_bytes() == unified.read_bytes()

    ref_asso = read_csv(REFERENCE / "asso_df.1trajectory.nknots4.csv")
    pkg_asso = read_csv(PACKAGED / "packaged_asso_df.1trajectory.nknots4.csv")
    summary["packaged_vs_historical_association"] = numeric_summary(ref_asso, pkg_asso)

    ref_sig = read_csv(REFERENCE / "sig_df.1trajectory.nknots4.csv")
    pkg_sig = read_csv(PACKAGED / "packaged_sig_df.1trajectory.nknots4.csv")
    summary["packaged_vs_historical_significant"] = numeric_summary(ref_sig, pkg_sig)
    summary["packaged_dynamic_call_same_genes"] = ref_sig.index.tolist() == pkg_sig.index.tolist()

    if (UNIFIED / "unified_asso_df.1trajectory.nknots4.csv").exists():
        uni_asso = read_csv(UNIFIED / "unified_asso_df.1trajectory.nknots4.csv")
        uni_sig = read_csv(UNIFIED / "unified_sig_df.1trajectory.nknots4.csv")
        summary["unified_vs_historical_association"] = numeric_summary(ref_asso, uni_asso)
        summary["unified_vs_historical_significant"] = numeric_summary(ref_sig, uni_sig)
        summary["unified_dynamic_call_same_genes"] = ref_sig.index.tolist() == uni_sig.index.tolist()

    for prefix, directory, asso_name, sig_name in [
        ("reference", REFERENCE, "asso_df.1trajectory.nknots4.csv", "sig_df.1trajectory.nknots4.csv"),
        ("packaged", PACKAGED, "packaged_asso_df.1trajectory.nknots4.csv", "packaged_sig_df.1trajectory.nknots4.csv"),
    ]:
        summary[f"{prefix}_file_rows"] = {
            "association": int(read_csv(directory / asso_name).shape[0]),
            "significant": int(read_csv(directory / sig_name).shape[0]),
        }
    if (UNIFIED / "unified_asso_df.1trajectory.nknots4.csv").exists():
        summary["unified_file_rows"] = {
            "association": int(read_csv(UNIFIED / "unified_asso_df.1trajectory.nknots4.csv").shape[0]),
            "significant": int(read_csv(UNIFIED / "unified_sig_df.1trajectory.nknots4.csv").shape[0]),
        }

    for prefix, directory, sce_name, gam_name in [
        ("reference", REFERENCE, "sce.1trajectory.nknots8.rds", "gamList.1trajectory.nknots8.rds"),
        ("packaged", PACKAGED, "packaged_sce.1trajectory.nknots8.rds", "packaged_gamList.1trajectory.nknots8.rds"),
    ]:
        summary[f"{prefix}_sce"] = rds_dims(directory / sce_name)
        summary[f"{prefix}_gamList"] = rds_dims(directory / gam_name)
    summary["packaged_vs_historical_sce_genes"] = rds_name_compare(
        REFERENCE / "sce.1trajectory.nknots8.rds",
        PACKAGED / "packaged_sce.1trajectory.nknots8.rds",
    )
    summary["packaged_vs_historical_gamList_genes"] = rds_name_compare(
        REFERENCE / "gamList.1trajectory.nknots8.rds",
        PACKAGED / "packaged_gamList.1trajectory.nknots8.rds",
    )
    summary["packaged_vs_historical_predictSmooth"] = predict_smooth_compare(
        REFERENCE / "gamList.1trajectory.nknots8.rds",
        PACKAGED / "packaged_gamList.1trajectory.nknots8.rds",
        n_points=9,
    )
    if (UNIFIED / "unified_sce.1trajectory.nknots8.rds").exists():
        summary["unified_sce"] = rds_dims(UNIFIED / "unified_sce.1trajectory.nknots8.rds")
        summary["unified_gamList"] = rds_dims(UNIFIED / "unified_gamList.1trajectory.nknots8.rds")
        summary["unified_vs_historical_sce_genes"] = rds_name_compare(
            REFERENCE / "sce.1trajectory.nknots8.rds",
            UNIFIED / "unified_sce.1trajectory.nknots8.rds",
        )
        summary["unified_vs_historical_gamList_genes"] = rds_name_compare(
            REFERENCE / "gamList.1trajectory.nknots8.rds",
            UNIFIED / "unified_gamList.1trajectory.nknots8.rds",
        )
        summary["unified_vs_historical_predictSmooth"] = predict_smooth_compare(
            REFERENCE / "gamList.1trajectory.nknots8.rds",
            UNIFIED / "unified_gamList.1trajectory.nknots8.rds",
            n_points=9,
        )

    out = RUN / "results" / "gene_dynamics_comparison_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
