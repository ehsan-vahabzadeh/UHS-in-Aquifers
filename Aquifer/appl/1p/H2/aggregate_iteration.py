import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--iter-id", required=True)
    args = parser.parse_args()

    iter_dir = Path("cases") / args.iter_id
    summaries = []

    for summary_path in iter_dir.rglob("case_summary.json"):
        with open(summary_path) as f:
            summaries.append(json.load(f))

    summaries = sorted(summaries, key=lambda x: x["case_id"])

    out_csv = Path("results") / f"{args.iter_id}_summary.csv"
    out_json = Path("results") / f"{args.iter_id}_summary.json"

    if summaries:
        fieldnames = sorted({
            key for s in summaries for key in s.keys()
            if key != "inputs"
        })

        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for s in summaries:
                flat = {k: s.get(k, "") for k in fieldnames}
                writer.writerow(flat)

        with open(out_json, "w") as f:
            json.dump(summaries, f, indent=2)

    print(f"Aggregated {len(summaries)} case summaries for {args.iter_id}")


if __name__ == "__main__":
    main()