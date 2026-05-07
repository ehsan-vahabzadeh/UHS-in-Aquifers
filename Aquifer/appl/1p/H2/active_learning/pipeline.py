"""
Pipeline for coupling DuMuX simulation summaries with a GP surrogate.

Components
----------
1. load_summary           -- read aggregate_iteration.py output CSV
2. wide_to_long           -- cycleXX_end_rf columns -> long rows with cycle_no
3. FeaturePrep            -- 8 continuous (normalized to [0,1]) + 4 one-hot cg_type
4. train_surrogate        -- BoTorch SingleTaskGP on 12-D input, scalar cumulative RF
5. score_candidates       -- predictive-variance acquisition, aggregated per case
6. diagnose_failures      -- simple comparison of failed vs successful inputs

The GP intentionally mirrors the kernel choices used in GP_ActiveLearning_v3.py
to stay consistent with the existing toy workflow.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# ------------------------------------------------------------------
# Schema -- continuous features and the canonical column-name aliases
# ------------------------------------------------------------------
CONTINUOUS_FEATURES = [
    # canonical name in our pipeline       -> list of acceptable source columns
    ("flow_rate_sm3_day",  ["FlowRateSm3Day",  "flow_rate_sm3_day"]),
    ("cycle_length_days",  ["CycleLengthDays", "cycle_length_days"]),
    ("cg_ratio",           ["CGRatio",         "cg_ratio"]),
    ("permeability_md",    ["PermeabilityMd",  "permeability_md"]),
    ("porosity",           ["Porosity",        "porosity"]),
    ("pressure_mpa",       ["PressureMPa",     "pressure_mpa"]),
    ("temperature_c",      ["TemperatureC",    "temperature_c"]),
]
CG_TYPES = ["CH4", "H2", "N2", "CO2"]
CG_TYPE_COLUMNS = ["CushionGasType", "cg_type"]
ONE_HOT_COLUMNS = [f"is_{g}" for g in CG_TYPES]
GP_FEATURE_COLUMNS = (
    [name for name, _ in CONTINUOUS_FEATURES] + ["cycle_no"] + ONE_HOT_COLUMNS
)
N_CONTINUOUS_GP = len(CONTINUOUS_FEATURES) + 1   # +1 for cycle_no
N_GP_FEATURES = N_CONTINUOUS_GP + len(ONE_HOT_COLUMNS)  # = 12

STATUS_COLUMNS = ["summary_status", "status"]


# ------------------------------------------------------------------
# 1. Loading & long-format conversion
# ------------------------------------------------------------------
def _pick_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def load_summary(csv_path: str | Path) -> pd.DataFrame:
    """Load an aggregate_iteration.py summary CSV with light cleanup."""
    df = pd.read_csv(csv_path, low_memory=False)
    return df


def _status_of(row: pd.Series) -> str:
    for c in STATUS_COLUMNS:
        if c in row and pd.notna(row[c]):
            return str(row[c]).strip().lower()
    return ""


def wide_to_long(
    df: pd.DataFrame,
    max_cycles: int = 10,
    target_prefix: str = "cycle",
    require_success: bool = True,
) -> pd.DataFrame:
    """
    Transform wide cycleXX_end_rf columns into long-format training rows.

    Returned DataFrame has the canonical input columns + ``cycle_no`` + ``cumulative_rf``.
    Rows with missing/blank target are skipped.
    """
    cont_cols = {}
    for canonical, aliases in CONTINUOUS_FEATURES:
        src = _pick_column(df, aliases)
        if src is None:
            raise ValueError(f"None of {aliases} are present in summary CSV")
        cont_cols[canonical] = src

    cg_src = _pick_column(df, CG_TYPE_COLUMNS)
    if cg_src is None:
        raise ValueError(f"None of {CG_TYPE_COLUMNS} are present in summary CSV")

    case_id_col = _pick_column(df, ["case_id", "summary_case_id"])

    rows = []
    for _, r in df.iterrows():
        if require_success and _status_of(r) != "success":
            continue

        base = {canonical: pd.to_numeric(r.get(src), errors="coerce")
                for canonical, src in cont_cols.items()}
        if any(pd.isna(v) for v in base.values()):
            continue
        cg = str(r.get(cg_src, "")).strip()
        if cg not in CG_TYPES:
            continue
        case_id = r.get(case_id_col) if case_id_col else None

        for cyc in range(1, max_cycles + 1):
            col = f"{target_prefix}{cyc:02d}_end_rf"
            if col not in df.columns:
                continue
            val = pd.to_numeric(r.get(col), errors="coerce")
            if pd.isna(val) or val == "":
                continue
            row = dict(base)
            row["cg_type"] = cg
            row["cycle_no"] = cyc
            row["cumulative_rf"] = float(val)
            if case_id is not None:
                row["case_id"] = case_id
            rows.append(row)

    return pd.DataFrame(rows)


# ------------------------------------------------------------------
# 2. Feature preparation -- continuous min/max normalization + one-hot
# ------------------------------------------------------------------
@dataclass
class FeaturePrep:
    """Fit min/max for continuous features and apply identical transform later."""
    bounds: dict = field(default_factory=dict)   # {col: (lo, hi)}
    max_cycles: int = 10
    feature_columns: list = field(default_factory=lambda: list(GP_FEATURE_COLUMNS))

    def fit(self, df: pd.DataFrame, candidate_df: pd.DataFrame | None = None,
            max_cycles: int | None = None) -> "FeaturePrep":
        if max_cycles is not None:
            self.max_cycles = max_cycles
        for canonical, _ in CONTINUOUS_FEATURES:
            vals = [pd.to_numeric(df[canonical], errors="coerce").dropna()]
            if candidate_df is not None and canonical in candidate_df.columns:
                vals.append(pd.to_numeric(candidate_df[canonical], errors="coerce").dropna())
            cat = pd.concat(vals, ignore_index=True)
            if cat.empty:
                raise ValueError(f"No usable values for {canonical} during FeaturePrep.fit")
            lo, hi = float(cat.min()), float(cat.max())
            if hi <= lo:
                hi = lo + 1.0
            self.bounds[canonical] = (lo, hi)
        # cycle_no bounded by [1, max_cycles]
        self.bounds["cycle_no"] = (1.0, float(max(self.max_cycles, 2)))
        return self

    def transform(self, df: pd.DataFrame) -> np.ndarray:
        cols = []
        for canonical, _ in CONTINUOUS_FEATURES:
            lo, hi = self.bounds[canonical]
            v = pd.to_numeric(df[canonical], errors="coerce").to_numpy(dtype=float)
            cols.append(np.clip((v - lo) / (hi - lo), 0.0, 1.0))
        lo, hi = self.bounds["cycle_no"]
        cyc = pd.to_numeric(df["cycle_no"], errors="coerce").to_numpy(dtype=float)
        cols.append(np.clip((cyc - lo) / (hi - lo), 0.0, 1.0))
        cg = df["cg_type"].astype(str).str.strip()
        for g in CG_TYPES:
            cols.append((cg == g).to_numpy(dtype=float))
        return np.column_stack(cols)


# ------------------------------------------------------------------
# 3. GP surrogate -- thin wrapper around BoTorch SingleTaskGP
# ------------------------------------------------------------------
def _import_gp_stack():
    import torch
    import gpytorch
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.input import Normalize
    from botorch.models.transforms.outcome import Standardize
    from gpytorch.mlls import ExactMarginalLogLikelihood
    from gpytorch.kernels import ScaleKernel, MaternKernel
    return (torch, gpytorch, fit_gpytorch_mll, SingleTaskGP,
            Normalize, Standardize, ExactMarginalLogLikelihood,
            ScaleKernel, MaternKernel)


def train_surrogate(X: np.ndarray, y: np.ndarray):
    """Train a SingleTaskGP. Inputs are already in [0,1]^12; targets are scalars."""
    (torch, _gp, fit_gpytorch_mll, SingleTaskGP,
     Normalize, Standardize, ExactMarginalLogLikelihood,
     ScaleKernel, MaternKernel) = _import_gp_stack()

    if X.ndim != 2 or X.shape[1] != N_GP_FEATURES:
        raise ValueError(f"Expected X with shape (n,{N_GP_FEATURES}), got {X.shape}")

    dtype = torch.double
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    train_X = torch.as_tensor(X, dtype=dtype, device=device)
    train_Y = torch.as_tensor(y, dtype=dtype, device=device).reshape(-1, 1)

    bounds = torch.stack([torch.zeros(N_GP_FEATURES, dtype=dtype, device=device),
                          torch.ones(N_GP_FEATURES, dtype=dtype, device=device)])
    covar = ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=N_GP_FEATURES))

    model = SingleTaskGP(
        train_X=train_X, train_Y=train_Y,
        covar_module=covar,
        input_transform=Normalize(d=N_GP_FEATURES, bounds=bounds),
        outcome_transform=Standardize(m=1),
    )
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    model.eval(); model.likelihood.eval()
    return model


def predict_mean_std(model, X: np.ndarray):
    import torch, gpytorch
    dtype = torch.double
    device = next(model.parameters()).device
    Xt = torch.as_tensor(X, dtype=dtype, device=device)
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        post = model.posterior(Xt)
        mean = post.mean.squeeze(-1).cpu().numpy()
        std = post.variance.clamp_min(0.0).sqrt().squeeze(-1).cpu().numpy()
    return mean, std


# ------------------------------------------------------------------
# 4. Candidate handling
# ------------------------------------------------------------------
def _normalize_candidate_df(pool_df: pd.DataFrame) -> pd.DataFrame:
    """Bring pool CSV columns to canonical names. Adds pool_row_index if missing."""
    df = pool_df.copy()
    for canonical, aliases in CONTINUOUS_FEATURES:
        if canonical not in df.columns:
            for a in aliases:
                if a in df.columns:
                    df[canonical] = pd.to_numeric(df[a], errors="coerce")
                    break
    if "cg_type" not in df.columns:
        for a in CG_TYPE_COLUMNS:
            if a in df.columns:
                df["cg_type"] = df[a].astype(str).str.strip()
                break
    if "pool_row_index" not in df.columns:
        df = df.reset_index(drop=True)
        df["pool_row_index"] = df.index
    return df


def exclude_already_run(
    candidates_df: pd.DataFrame, summary_df: pd.DataFrame
) -> tuple[pd.DataFrame, int]:
    """Drop candidates whose pool_row_index already appears in the summary."""
    if "pool_row_index" not in summary_df.columns:
        return candidates_df, 0
    done = set(pd.to_numeric(summary_df["pool_row_index"], errors="coerce").dropna().astype(int))
    keep = ~candidates_df["pool_row_index"].astype(int).isin(done)
    return candidates_df[keep].reset_index(drop=True), int((~keep).sum())


def expand_candidates_with_cycles(
    cand_df: pd.DataFrame, max_cycles: int
) -> pd.DataFrame:
    """For each physical case create one row per cycle_no in 1..max_cycles."""
    rows = []
    for _, r in cand_df.iterrows():
        for cyc in range(1, max_cycles + 1):
            row = r.to_dict()
            row["cycle_no"] = cyc
            rows.append(row)
    return pd.DataFrame(rows)


def score_candidates(
    model, prep: FeaturePrep, cand_df: pd.DataFrame,
    max_cycles: int, aggregate: str = "mean",
) -> pd.DataFrame:
    """
    Score every candidate by acquisition value.

    Returns a per-physical-case DataFrame with score (and per-cycle stats),
    sorted descending by score.
    """
    expanded = expand_candidates_with_cycles(cand_df, max_cycles)
    X = prep.transform(expanded)
    mean, std = predict_mean_std(model, X)
    expanded = expanded.copy()
    expanded["pred_mean"] = mean
    expanded["pred_std"] = std

    grp = expanded.groupby("pool_row_index", sort=False)
    if aggregate == "max":
        scores = grp["pred_std"].max()
    elif aggregate == "sum":
        scores = grp["pred_std"].sum()
    else:
        scores = grp["pred_std"].mean()

    summary = (cand_df.set_index("pool_row_index").loc[scores.index].reset_index())
    summary["score"] = scores.values
    summary["pred_std_mean"] = grp["pred_std"].mean().values
    summary["pred_std_max"] = grp["pred_std"].max().values
    summary["pred_mean_mean"] = grp["pred_mean"].mean().values
    return summary.sort_values("score", ascending=False).reset_index(drop=True)


# ------------------------------------------------------------------
# 5. Failure diagnostics
# ------------------------------------------------------------------
def diagnose_failures(df: pd.DataFrame) -> dict:
    """
    Compare per-feature ranges of failed vs successful simulations.

    Returns a small dict suitable for printing.
    """
    status = df.apply(_status_of, axis=1)
    succ = df[status == "success"]
    fail = df[(status != "") & (status != "success")]
    out = {"n_total": int(len(df)),
           "n_success": int(len(succ)),
           "n_failed": int(len(fail)),
           "per_feature": {}}
    if fail.empty or succ.empty:
        return out
    for canonical, aliases in CONTINUOUS_FEATURES:
        src = _pick_column(df, aliases)
        if src is None:
            continue
        s = pd.to_numeric(succ[src], errors="coerce").dropna()
        f = pd.to_numeric(fail[src], errors="coerce").dropna()
        if s.empty or f.empty:
            continue
        out["per_feature"][canonical] = {
            "succ_min": float(s.min()), "succ_max": float(s.max()),
            "fail_min": float(f.min()), "fail_max": float(f.max()),
            "fail_mean": float(f.mean()), "succ_mean": float(s.mean()),
        }
    cg_src = _pick_column(df, CG_TYPE_COLUMNS)
    if cg_src is not None:
        out["fail_cg_counts"] = fail[cg_src].astype(str).value_counts().to_dict()
        out["succ_cg_counts"] = succ[cg_src].astype(str).value_counts().to_dict()
    return out
