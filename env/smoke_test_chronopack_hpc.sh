#!/usr/bin/env bash
set -euo pipefail

repo="${CHRONOPACK_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
tmp_root="${TMPDIR:-/tmp}"
mkdir -p "${tmp_root}/chronopack-smoke"

echo "== Python package import from outside repository =="
(
  cd "${tmp_root}/chronopack-smoke"
  python - <<'PY'
from importlib import resources

import anndata
import gseapy
import matplotlib
import networkx
import numpy
import omnipath
import pandas
import scanpy
import scipy
import seaborn
import sklearn
import statsmodels
import wot

import chronopack as cp
from chronopack._backends.r import available_r_resources, read_r_resource_text

assert cp.io is not None
assert cp.tl is not None
assert cp.pl is not None
assert "tradeseq_helpers.R" in available_r_resources()
assert "fit_tradeseq_sce <- function" in read_r_resource_text("tradeseq_helpers.R")
assert resources.files("chronopack._backends.r").joinpath("tradeseq_helpers.R").is_file()

print("chronopack", cp.__file__)
print("resources", ", ".join(available_r_resources()))
PY
)

echo "== WOT toy transport smoke test =="
python - <<'PY'
import anndata as ad
import numpy as np
import pandas as pd
import wot

adata = ad.AnnData(
    X=np.array([[1.0, 0.0], [0.0, 1.0], [0.9, 0.1], [0.1, 0.9]]),
    obs=pd.DataFrame({"day": [0.0, 0.0, 1.0, 1.0]}, index=["a", "b", "c", "d"]),
)
model = wot.ot.OTModel(adata, day_field="day")
tm = model.compute_transport_map(0.0, 1.0)
assert tm.X.shape == (2, 2)
assert np.isfinite(tm.X).all()
print("transport_map_shape", tm.X.shape)
PY

echo "== R backend packages =="
Rscript - <<'RS'
suppressPackageStartupMessages({
  library(tradeSeq)
  library(SingleCellExperiment)
  library(fdrtool)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(matrixStats)
  library(RColorBrewer)
  library(patchwork)
  library(circlize)
  library(cluster)
  library(ggrepel)
  library(shadowtext)
  library(purrr)
  library(scales)
  library(slingshot)
  library(BiocParallel)
  library(SummarizedExperiment)
})
cat("R.version:", R.version.string, "\n")
RS

echo "== rpy2 Python -> R -> Python round trip =="
python - <<'PY'
from rpy2.robjects import FloatVector, r

values = FloatVector([1.0, 2.0, 3.0])
mean_value = float(r["mean"](values)[0])
assert mean_value == 2.0
print("rpy2_mean", mean_value)
PY

echo "== pytest =="
cd "${repo}"
python -m pytest
