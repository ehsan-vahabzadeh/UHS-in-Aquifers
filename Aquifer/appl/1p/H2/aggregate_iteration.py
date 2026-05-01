import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--iter-id", required=True)
    args = parser.parse_args()

    with open(args.manifest, newline="") as f:
        manifest_rows = list(csv.DictReader(f))

    iter_dir = Path("cases") / args.iter_id
    summaries = []

    for summary_path in iter_dir.rglob("case_summary.json"):
        with open(summary_path) as f:
            summaries.append(json.load(f))

    summaries = sorted(summaries, key=lambda x: x["case_id"])

    out_csv = Path("results") / f"{args.iter_id}_summary.csv"
    out_json = Path("results") / f"{args.iter_id}_summary.json"
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = sorted({
        key for s in summaries for key in s.keys()
        if key != "inputs"
    })
    if not fieldnames:
        fieldnames = ["case_id", "chunk_id", "name", "status", "case_dir", "exception", "returncode"]

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in summaries:
            flat = {k: s.get(k, "") for k in fieldnames}
            writer.writerow(flat)

    with open(out_json, "w") as f:
        json.dump(summaries, f, indent=2)

    expected = len(manifest_rows)
    found = len(summaries)
    print(f"Manifest rows expected: {expected}")
    print(f"Case summaries found: {found}")
    print(f"Wrote CSV summary: {out_csv}")
    print(f"Wrote JSON summary: {out_json}")

    if found != expected:
        missing = sorted(
            set(row["case_id"] for row in manifest_rows)
            - set(str(summary["case_id"]) for summary in summaries)
        )
        print(f"WARNING: missing case summaries for case_id values: {', '.join(missing)}")

    print(f"Aggregated {found} case summaries for {args.iter_id}")


if __name__ == "__main__":
    main()
