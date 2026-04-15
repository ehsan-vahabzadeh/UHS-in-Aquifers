#!/usr/bin/env python3
"""Generate an initial two-phase Aquifer simulation campaign and input files.

Design philosophy:
1) Cluster real reservoirs in physical-property space
2) Allocate case budget across clusters with minimum representation + heterogeneity weighting
3) Sample operational variables with space-filling design (LHS or Sobol)
4) Balance cushion gas assignment across cases
5) Emit manifest CSV/JSON and one simulator input file per case

This script is intentionally extensible for later active-learning rounds (append mode).
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from scipy.stats import qmc
except Exception:  # pragma: no cover - optional dependency at runtime
    qmc = None

try:
    from sklearn.cluster import KMeans
except Exception:  # pragma: no cover - optional dependency at runtime
    KMeans = None


LOGGER = logging.getLogger("aquifer_campaign_generator")
GAS_TYPES = ("CO2", "N2", "CH4", "H2")


@dataclass(frozen=True)
class ReservoirRecord:
    reservoir_id: str
    porosity: float
    permeability: float  # mD expected
    pressure: float  # MPa expected
    temperature: float  # Celsius expected
    cluster_id: int


@dataclass(frozen=True)
class CampaignCase:
    case_id: str
    reservoir_id: str
    cluster_id: int
    porosity: float
    permeability: float
    pressure: float
    temperature: float
    cushion_gas_type: str
    flow_rate: float
    cushion_gas_ratio: float
    cycle_length: float
    injection_duration_dev: float
    injection_duration_op: float
    extraction_duration_op: float
    t_end_days: float
    injection_rate_dev: float
    injection_rate_op: float
    production_rate: float
    reference_porosity: float
    reference_permeability_m2: float
    pressure_top_pa: float
    initial_temperature_k: float
    input_file: str
    notes: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reservoir-csv", type=Path, required=True)
    parser.add_argument("--template-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-cases", type=int, default=300)
    parser.add_argument("--n-clusters", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--design-method", choices=("lhs", "sobol"), default="lhs")
    parser.add_argument("--cases-per-reservoir-min", type=int, default=1)
    parser.add_argument("--cases-per-reservoir-max", type=int, default=20)
    parser.add_argument("--per-gas-balance", action="store_true", help="Balance gas types globally")
    parser.add_argument("--append", action="store_true", help="Append cases to existing campaign_manifest.csv")
    parser.add_argument("--manifest-name", default="campaign_manifest.csv")
    parser.add_argument("--cluster-weight-size", type=float, default=0.6)
    parser.add_argument("--cluster-weight-spread", type=float, default=0.4)
    parser.add_argument("--flow-min", type=float, default=100000.0)
    parser.add_argument("--flow-max", type=float, default=1500000.0)
    parser.add_argument("--ratio-min", type=float, default=0.0)
    parser.add_argument("--ratio-max", type=float, default=3.0)
    parser.add_argument("--cycle-min", type=float, default=14.0)
    parser.add_argument("--cycle-max", type=float, default=365.0)
    parser.add_argument("--n-operation-cycles", type=float, default=10.0)
    parser.add_argument("--idle-op-days", type=float, default=0.0)
    parser.add_argument("--well-height", type=float, default=10.0)
    parser.add_argument("--well-radius", type=float, default=0.2)
    parser.add_argument(
        "--sm3-to-mol-factor",
        type=float,
        default=0.041e3,
        help="Conversion factor from Sm3/day to mol/day at standard conditions",
    )
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    return parser.parse_args()


def standardize(arr: np.ndarray) -> np.ndarray:
    mean = np.mean(arr, axis=0)
    std = np.std(arr, axis=0)
    std[std == 0.0] = 1.0
    return (arr - mean) / std


def load_reservoirs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"reservoir_id", "porosity", "permeability", "pressure", "temperature"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required reservoir columns: {sorted(missing)}")

    for col in ["porosity", "permeability", "pressure", "temperature"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["reservoir_id", "porosity", "permeability", "pressure", "temperature"]).copy()
    if (df["permeability"] <= 0).any():
        raise ValueError("Permeability must be strictly positive for log transform.")
    df["reservoir_id"] = df["reservoir_id"].astype(str)
    return df.reset_index(drop=True)


def cluster_reservoirs(df: pd.DataFrame, n_clusters: int, seed: int) -> Tuple[pd.DataFrame, np.ndarray]:
    x = np.column_stack(
        [
            df["porosity"].to_numpy(),
            np.log(df["permeability"].to_numpy()),
            df["pressure"].to_numpy(),
            df["temperature"].to_numpy(),
        ]
    )
    x_std = standardize(x)

    k = max(1, min(n_clusters, len(df)))
    if KMeans is not None:
        model = KMeans(n_clusters=k, random_state=seed, n_init="auto")
        labels = model.fit_predict(x_std)
        centers = model.cluster_centers_
    else:
        LOGGER.warning("scikit-learn unavailable; using deterministic fallback clustering.")
        sort_idx = np.argsort(x_std[:, 1])
        labels = np.zeros(len(df), dtype=int)
        for i, idx in enumerate(sort_idx):
            labels[idx] = int(i * k / len(df))
        centers = np.vstack([x_std[labels == cid].mean(axis=0) for cid in range(k)])

    out = df.copy()
    out["cluster_id"] = labels
    return out, centers


def cluster_spread(x_std: np.ndarray, labels: np.ndarray, centers: np.ndarray) -> Dict[int, float]:
    spread: Dict[int, float] = {}
    for cid in range(centers.shape[0]):
        cluster_points = x_std[labels == cid]
        if len(cluster_points) == 0:
            spread[cid] = 0.0
            continue
        d = np.linalg.norm(cluster_points - centers[cid], axis=1)
        spread[cid] = float(np.mean(d))
    return spread


def allocate_budget(
    df: pd.DataFrame,
    total_cases: int,
    size_weight: float,
    spread_weight: float,
    cpr_min: int,
    cpr_max: int,
) -> Dict[int, int]:
    counts = df.groupby("cluster_id").size().to_dict()
    cids = sorted(counts)
    min_per_cluster = 1
    if total_cases < len(cids) * min_per_cluster:
        raise ValueError("n-cases too small to give each cluster minimum representation.")

    x_std = standardize(
        np.column_stack(
            [
                df["porosity"].to_numpy(),
                np.log(df["permeability"].to_numpy()),
                df["pressure"].to_numpy(),
                df["temperature"].to_numpy(),
            ]
        )
    )
    labels = df["cluster_id"].to_numpy()
    centers = np.vstack([x_std[labels == cid].mean(axis=0) for cid in cids])
    spreads = cluster_spread(x_std, labels, centers)

    size_arr = np.array([counts[c] for c in cids], dtype=float)
    spread_arr = np.array([spreads[c] for c in cids], dtype=float)
    if size_arr.sum() == 0:
        size_arr[:] = 1.0
    if spread_arr.sum() == 0:
        spread_arr[:] = 1.0

    size_norm = size_arr / size_arr.sum()
    spread_norm = spread_arr / spread_arr.sum()
    weights = size_weight * size_norm + spread_weight * spread_norm
    weights = weights / weights.sum()

    remaining = total_cases - len(cids) * min_per_cluster
    raw = remaining * weights
    add = np.floor(raw).astype(int)
    leftover = int(remaining - add.sum())
    order = np.argsort(-(raw - add))
    for i in range(leftover):
        add[order[i % len(order)]] += 1

    alloc = {cid: int(min_per_cluster + add[i]) for i, cid in enumerate(cids)}

    # Cap by representative reservoir max if requested.
    for cid in cids:
        n_res = int(counts[cid])
        alloc[cid] = min(alloc[cid], max(cpr_min, n_res * cpr_max))

    # Re-distribute if capping removed cases.
    deficit = total_cases - sum(alloc.values())
    if deficit > 0:
        for cid in cids:
            if deficit == 0:
                break
            alloc[cid] += 1
            deficit -= 1

    return alloc


def select_representatives(df: pd.DataFrame, seed: int) -> Dict[int, List[ReservoirRecord]]:
    rng = np.random.default_rng(seed)
    reps: Dict[int, List[ReservoirRecord]] = {}

    for cid, gdf in df.groupby("cluster_id"):
        g = gdf.copy()
        g["_dist"] = np.sqrt(
            (standardize(g[["porosity"]].to_numpy())[:, 0] ** 2)
            + (standardize(np.log(g[["permeability"]].to_numpy()))[:, 0] ** 2)
            + (standardize(g[["pressure"]].to_numpy())[:, 0] ** 2)
            + (standardize(g[["temperature"]].to_numpy())[:, 0] ** 2)
        )
        g = g.sort_values("_dist")
        ordered = list(g.to_dict("records"))
        if len(ordered) > 1:
            # small shuffle among top candidates to avoid deterministic overuse
            top_n = min(5, len(ordered))
            idx = np.arange(top_n)
            rng.shuffle(idx)
            top = [ordered[i] for i in idx]
            ordered = top + ordered[top_n:]

        reps[cid] = [
            ReservoirRecord(
                reservoir_id=str(r["reservoir_id"]),
                porosity=float(r["porosity"]),
                permeability=float(r["permeability"]),
                pressure=float(r["pressure"]),
                temperature=float(r["temperature"]),
                cluster_id=int(cid),
            )
            for r in ordered
        ]

    return reps


def qmc_unit_samples(n: int, dim: int, method: str, seed: int) -> np.ndarray:
    if qmc is None:
        rng = np.random.default_rng(seed)
        return rng.random((n, dim))

    if method == "sobol":
        sampler = qmc.Sobol(d=dim, scramble=True, seed=seed)
        m = math.ceil(math.log2(max(n, 2)))
        u = sampler.random_base2(m=m)
        return u[:n]

    sampler = qmc.LatinHypercube(d=dim, seed=seed)
    return sampler.random(n=n)


def balanced_gases(n: int, seed: int) -> List[str]:
    base = [GAS_TYPES[i % len(GAS_TYPES)] for i in range(n)]
    rng = np.random.default_rng(seed)
    rng.shuffle(base)
    return base


def flow_to_molar_rate(flow_sm3_day: float, sm3_to_mol_factor: float, well_height: float, well_radius: float) -> float:
    """Approximate boundary molar flux [mol m-2 s-1] used by existing Aquifer scripts.

    Assumption retained from existing generator style:
    - convert Sm3/day -> mol/day by factor
    - divide by 86400
    - divide by circumferential flow area proxy (2*pi*r*H)
    """
    area = max(1e-9, 2.0 * math.pi * well_radius * well_height)
    return (flow_sm3_day * sm3_to_mol_factor) / 86400.0 / area


def sanitize(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


def build_cases(
    reps: Dict[int, List[ReservoirRecord]],
    allocation: Dict[int, int],
    args: argparse.Namespace,
) -> List[CampaignCase]:
    rng = np.random.default_rng(args.seed)
    cases: List[CampaignCase] = []

    # Balanced gas plan at campaign level (optional)
    if args.per_gas_balance:
        gas_pool = balanced_gases(sum(allocation.values()), args.seed)
    else:
        gas_pool = [GAS_TYPES[i] for i in rng.integers(0, len(GAS_TYPES), size=sum(allocation.values()))]

    gas_ptr = 0
    for cid in sorted(allocation):
        n_cases = allocation[cid]
        if n_cases <= 0:
            continue
        u = qmc_unit_samples(n_cases, 3, args.design_method, args.seed + cid)
        cluster_reps = reps[cid]

        for i in range(n_cases):
            rec = cluster_reps[i % len(cluster_reps)]
            flow = args.flow_min + u[i, 0] * (args.flow_max - args.flow_min)
            ratio = args.ratio_min + u[i, 1] * (args.ratio_max - args.ratio_min)
            cycle = args.cycle_min + u[i, 2] * (args.cycle_max - args.cycle_min)
            gas = gas_pool[gas_ptr]
            gas_ptr += 1

            inj_op = cycle / 2.0
            ext_op = cycle / 2.0
            inj_dev = ratio * inj_op
            t_end_days = args.n_operation_cycles * cycle + inj_dev

            molar_flux = flow_to_molar_rate(
                flow_sm3_day=flow,
                sm3_to_mol_factor=args.sm3_to_mol_factor,
                well_height=args.well_height,
                well_radius=args.well_radius,
            )

            case_id = f"c{len(cases)+1:05d}_{sanitize(rec.reservoir_id)}_{gas}"
            case = CampaignCase(
                case_id=case_id,
                reservoir_id=rec.reservoir_id,
                cluster_id=rec.cluster_id,
                porosity=rec.porosity,
                permeability=rec.permeability,
                pressure=rec.pressure,
                temperature=rec.temperature,
                cushion_gas_type=gas,
                flow_rate=flow,
                cushion_gas_ratio=ratio,
                cycle_length=cycle,
                injection_duration_dev=inj_dev,
                injection_duration_op=inj_op,
                extraction_duration_op=ext_op,
                t_end_days=t_end_days,
                injection_rate_dev=-molar_flux,
                injection_rate_op=-molar_flux,
                production_rate=molar_flux,
                reference_porosity=rec.porosity,
                reference_permeability_m2=rec.permeability * 9.86923e-16,
                pressure_top_pa=rec.pressure * 1e6,
                initial_temperature_k=rec.temperature + 273.15,
                input_file=str(Path("inputs") / f"{case_id}.input"),
            )
            cases.append(case)

    return cases


def update_params_text(
    template_text: str,
    section_key_values: Dict[str, Dict[str, str]],
    add_missing_keys: bool = True,
) -> str:
    lines = template_text.splitlines()
    current_section: Optional[str] = None
    seen: Dict[Tuple[str, str], bool] = {}

    section_pat = re.compile(r"^\s*\[([^\]]+)\]\s*$")
    kv_pat = re.compile(r"^(\s*)([A-Za-z0-9_.]+)(\s*=\s*)([^#]*?)(\s*(#.*)?)$")

    out_lines: List[str] = []
    for line in lines:
        sec_match = section_pat.match(line)
        if sec_match:
            current_section = sec_match.group(1)
            out_lines.append(line)
            continue

        m = kv_pat.match(line)
        if m and current_section is not None:
            key = m.group(2)
            if current_section in section_key_values and key in section_key_values[current_section]:
                new_val = section_key_values[current_section][key]
                out_lines.append(f"{m.group(1)}{key}{m.group(3)}{new_val}{m.group(5)}")
                seen[(current_section, key)] = True
            else:
                out_lines.append(line)
        else:
            out_lines.append(line)

    if add_missing_keys:
        # append any missing keys under existing section or at end if section absent
        for section, mapping in section_key_values.items():
            section_exists = any(section_pat.match(ln) and section_pat.match(ln).group(1) == section for ln in out_lines)
            if not section_exists:
                out_lines.append("")
                out_lines.append(f"[{section}]")

            insert_idx = len(out_lines)
            for i, ln in enumerate(out_lines):
                sm = section_pat.match(ln)
                if sm and sm.group(1) == section:
                    insert_idx = i + 1
                    while insert_idx < len(out_lines) and not section_pat.match(out_lines[insert_idx]):
                        insert_idx += 1
                    break

            new_items = []
            for key, val in mapping.items():
                if (section, key) not in seen:
                    new_items.append(f"{key} = {val}")
            if new_items:
                out_lines[insert_idx:insert_idx] = new_items

    return "\n".join(out_lines) + "\n"


def case_to_param_map(case: CampaignCase) -> Dict[str, Dict[str, str]]:
    # Keep mapping conservative and aligned with observed Aquifer *.input style.
    mapping: Dict[str, Dict[str, str]] = {
        "Problem": {
            "Name": case.case_id,
            "InitialTemperature": f"{case.initial_temperature_k:.6f}",
        },
        "BoundaryConditions": {
            "Pressure_TOP": f"{case.pressure_top_pa:.6f}",
            "InjectionDurationDev": f"{case.injection_duration_dev:.6f}",
            "InjectionDurationOp": f"{case.injection_duration_op:.6f}",
            "ExtractionDurationOp": f"{case.extraction_duration_op:.6f}",
            "IdleDurationOp": "0.0",
            "InjectionRateDev": f"{case.injection_rate_dev:.10f}",
            "InjectionRateOp": f"{case.injection_rate_op:.10f}",
            "ProductionRate": f"{case.production_rate:.10f}",
            "Well_Height": "10",
            "CushionGasType": case.cushion_gas_type,
        },
        "SpatialParams": {
            "ReferencePorosity": f"{case.reference_porosity:.6f}",
            "ReferencePermeability": f"{case.reference_permeability_m2:.10e}",
        },
        "TimeLoop": {
            "TEnd": f"{case.t_end_days:.6f}",
        },
    }

    # Compatibility with templates expecting HydrogenInjectionConcentration key.
    if case.cushion_gas_type == "H2":
        mapping["BoundaryConditions"]["HydrogenInjectionConcentration"] = "1"
    else:
        mapping["BoundaryConditions"]["HydrogenInjectionConcentration"] = "0"

    return mapping


def write_inputs(cases: Sequence[CampaignCase], template_file: Path, out_dir: Path) -> None:
    template_text = template_file.read_text(encoding="utf-8")
    inputs_dir = out_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    for case in cases:
        updated = update_params_text(template_text, case_to_param_map(case), add_missing_keys=True)
        (inputs_dir / f"{case.case_id}.input").write_text(updated, encoding="utf-8")


def diagnostics(cases: Sequence[CampaignCase], out_dir: Path) -> None:
    if not cases:
        LOGGER.warning("No cases generated.")
        return

    df = pd.DataFrame([asdict(c) for c in cases])
    by_cluster = df.groupby("cluster_id").size().to_dict()
    by_gas = df.groupby("cushion_gas_type").size().to_dict()

    LOGGER.info("Generated %d cases", len(cases))
    LOGGER.info("Cases per cluster: %s", by_cluster)
    LOGGER.info("Cases per gas: %s", by_gas)

    coverage = {
        "flow_rate": (float(df["flow_rate"].min()), float(df["flow_rate"].max())),
        "cushion_gas_ratio": (float(df["cushion_gas_ratio"].min()), float(df["cushion_gas_ratio"].max())),
        "cycle_length": (float(df["cycle_length"].min()), float(df["cycle_length"].max())),
        "porosity": (float(df["porosity"].min()), float(df["porosity"].max())),
        "permeability": (float(df["permeability"].min()), float(df["permeability"].max())),
        "pressure": (float(df["pressure"].min()), float(df["pressure"].max())),
        "temperature": (float(df["temperature"].min()), float(df["temperature"].max())),
    }

    report = {
        "n_cases": len(cases),
        "n_clusters": int(df["cluster_id"].nunique()),
        "n_reservoirs_used": int(df["reservoir_id"].nunique()),
        "cases_per_cluster": by_cluster,
        "cases_per_gas": by_gas,
        "coverage": coverage,
    }
    (out_dir / "campaign_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


def save_manifest(cases: Sequence[CampaignCase], out_dir: Path, manifest_name: str, append: bool) -> None:
    path = out_dir / manifest_name
    df = pd.DataFrame([asdict(c) for c in cases])

    if append and path.exists():
        prev = pd.read_csv(path)
        df = pd.concat([prev, df], ignore_index=True)
        df = df.drop_duplicates(subset=["case_id"], keep="last")

    df.to_csv(path, index=False)


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = load_reservoirs(args.reservoir_csv)
    LOGGER.info("Loaded %d reservoirs", len(df))

    clustered_df, _centers = cluster_reservoirs(df, args.n_clusters, args.seed)
    LOGGER.info("Created %d clusters", clustered_df['cluster_id'].nunique())

    reps = select_representatives(clustered_df, args.seed)
    alloc = allocate_budget(
        clustered_df,
        total_cases=args.n_cases,
        size_weight=args.cluster_weight_size,
        spread_weight=args.cluster_weight_spread,
        cpr_min=args.cases_per_reservoir_min,
        cpr_max=args.cases_per_reservoir_max,
    )

    cases = build_cases(reps, alloc, args)

    write_inputs(cases, args.template_file, args.output_dir)
    save_manifest(cases, args.output_dir, args.manifest_name, append=args.append)
    diagnostics(cases, args.output_dir)

    LOGGER.info("Wrote manifest: %s", args.output_dir / args.manifest_name)
    LOGGER.info("Wrote input files under: %s", args.output_dir / "inputs")


if __name__ == "__main__":
    main()
