#!/usr/bin/env python3
import argparse
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


PVTU_RE = re.compile(r"^s(?P<size>\d+)-(?P<simulation>.+)-(?P<step>\d+)\.pvtu$")
VTU_RE = re.compile(r"^s(?P<size>\d+)-p(?P<part>\d+)-(?P<simulation>.+)-(?P<step>\d+)\.vtu$")


@dataclass(frozen=True)
class Dataset:
    simulation: str
    timestep: str
    step: str
    num_parts: int
    source: Path | None
    pieces: tuple[Path, ...]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge parallel DuMuX VTU output into serial VTU files."
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        type=Path,
        help="Directory containing the .pvd/.pvtu/.vtu output files.",
    )
    parser.add_argument(
        "--num-parts",
        type=int,
        help="Expected MPI part count, usually the number passed to mpirun -np.",
    )
    parser.add_argument(
        "--simulation-name",
        help="Only merge output for this Problem.Name.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rewrite merged .vtu files that already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print what would be merged.",
    )
    return parser.parse_args()


def parse_pvtu_name(path):
    match = PVTU_RE.match(path.name)
    if not match:
        return None

    return {
        "num_parts": int(match.group("size")),
        "simulation": match.group("simulation"),
        "step": match.group("step"),
    }


def parse_vtu_name(path):
    match = VTU_RE.match(path.name)
    if not match:
        return None

    return {
        "num_parts": int(match.group("size")),
        "part": int(match.group("part")),
        "simulation": match.group("simulation"),
        "step": match.group("step"),
    }


def read_pvtu_pieces(pvtu_path):
    try:
        tree = ET.parse(pvtu_path)
    except ET.ParseError as exc:
        print(f"Warning: could not parse {pvtu_path}: {exc}")
        return ()

    pieces = []
    for piece in tree.findall(".//Piece"):
        source = piece.attrib.get("Source")
        if source:
            pieces.append((pvtu_path.parent / source).resolve())

    return tuple(pieces)


def datasets_from_pvd(pvd_path, expected_parts, simulation_name):
    try:
        tree = ET.parse(pvd_path)
    except ET.ParseError as exc:
        print(f"Warning: could not parse {pvd_path}: {exc}")
        return []

    datasets = []
    part_mismatch_counts = {}
    missing_source_count = 0
    missing_piece_dataset_count = 0
    missing_piece_file_count = 0
    for data_set in tree.findall(".//DataSet"):
        file_name = data_set.attrib.get("file")
        if not file_name:
            continue

        source = (pvd_path.parent / file_name).resolve()
        parsed = parse_pvtu_name(source)
        if not parsed:
            continue

        if simulation_name is not None and parsed["simulation"] != simulation_name:
            continue

        if expected_parts is not None and parsed["num_parts"] != expected_parts:
            part_mismatch_counts[parsed["num_parts"]] = (
                part_mismatch_counts.get(parsed["num_parts"], 0) + 1
            )
            continue

        if not source.exists():
            missing_source_count += 1
            continue

        pieces = read_pvtu_pieces(source)
        missing_pieces = [piece for piece in pieces if not piece.exists()]
        if missing_pieces:
            missing_piece_dataset_count += 1
            missing_piece_file_count += len(missing_pieces)
            continue

        if expected_parts is not None and len(pieces) != expected_parts:
            print(
                f"Skipping {source.name}: it references {len(pieces)} parts, "
                f"expected {expected_parts}."
            )
            continue

        datasets.append(
            Dataset(
                simulation=parsed["simulation"],
                timestep=data_set.attrib.get("timestep", parsed["step"]),
                step=parsed["step"],
                num_parts=parsed["num_parts"],
                source=source,
                pieces=pieces,
            )
        )

    for actual_parts, count in sorted(part_mismatch_counts.items()):
        print(
            f"Skipping {count} dataset(s) in {pvd_path.name}: written for "
            f"{actual_parts} parts, expected {expected_parts}."
        )

    if missing_source_count:
        print(
            f"Warning: skipping {missing_source_count} dataset(s) in {pvd_path.name}; "
            "the referenced .pvtu file is missing."
        )

    if missing_piece_dataset_count:
        print(
            f"Warning: skipping {missing_piece_dataset_count} dataset(s) in "
            f"{pvd_path.name}; {missing_piece_file_count} referenced piece file(s) "
            "are missing."
        )

    return datasets


def datasets_from_vtu_files(input_dir, expected_parts, simulation_name):
    grouped = {}
    for path in input_dir.glob("s*-p*-*.vtu"):
        parsed = parse_vtu_name(path)
        if not parsed:
            continue

        if simulation_name is not None and parsed["simulation"] != simulation_name:
            continue

        if expected_parts is not None and parsed["num_parts"] != expected_parts:
            continue

        key = (parsed["simulation"], parsed["step"], parsed["num_parts"])
        grouped.setdefault(key, {})[parsed["part"]] = path.resolve()

    datasets = []
    for (simulation, step, num_parts), pieces_by_part in sorted(grouped.items()):
        missing_parts = [
            part for part in range(num_parts) if part not in pieces_by_part
        ]
        if missing_parts:
            missing = ", ".join(f"p{part:04}" for part in missing_parts[:8])
            if len(missing_parts) > 8:
                missing += ", ..."
            print(
                f"Warning: skipping {simulation}-{step}; missing {len(missing_parts)} "
                f"of {num_parts} parts ({missing})."
            )
            continue

        pieces = tuple(pieces_by_part[part] for part in range(num_parts))
        datasets.append(
            Dataset(
                simulation=simulation,
                timestep=step,
                step=step,
                num_parts=num_parts,
                source=None,
                pieces=pieces,
            )
        )

    return datasets


def collect_datasets(input_dir, expected_parts, simulation_name):
    datasets = []
    for pvd_path in sorted(input_dir.glob("*.pvd")):
        if simulation_name is not None and pvd_path.stem != simulation_name:
            continue

        pvd_datasets = datasets_from_pvd(pvd_path, expected_parts, simulation_name)
        if pvd_datasets:
            print(f"Found {len(pvd_datasets)} timestep(s) in {pvd_path.name}.")
            datasets.extend(pvd_datasets)

    if datasets:
        return datasets

    print("No usable .pvd/.pvtu datasets found; scanning .vtu pieces directly.")
    return datasets_from_vtu_files(input_dir, expected_parts, simulation_name)


def merge_dataset(dataset, output_path):
    import vtk

    if dataset.source is not None:
        reader = vtk.vtkXMLPUnstructuredGridReader()
        reader.SetFileName(str(dataset.source))
        reader.Update()
        grid = reader.GetOutput()
    else:
        append_filter = vtk.vtkAppendFilter()
        for piece in dataset.pieces:
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(str(piece))
            reader.Update()
            append_filter.AddInputData(reader.GetOutput())
        append_filter.Update()
        grid = append_filter.GetOutput()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputData(grid)
    if writer.Write() != 1:
        raise RuntimeError(f"Could not write merged output {output_path}")


def write_pvd(output_dir, simulation, datasets):
    pvd_path = output_dir / f"{simulation}.pvd"
    with pvd_path.open("w", encoding="utf-8") as pvd:
        pvd.write('<?xml version="1.0"?>\n')
        pvd.write(
            '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
        )
        pvd.write("<Collection>\n")
        for dataset in datasets:
            file_name = merged_file_name(dataset)
            pvd.write(
                f'<DataSet timestep="{dataset.timestep}" group="" part="0" '
                f'name="" file="{file_name}"/>\n'
            )
        pvd.write("</Collection>\n")
        pvd.write("</VTKFile>\n")


def merged_file_name(dataset):
    return f"s{dataset.num_parts:04}-{dataset.simulation}-{dataset.step}.vtu"


def delete_parallel_files(input_dir, datasets):
    files_to_delete = set()
    for dataset in datasets:
        if dataset.source is not None:
            files_to_delete.add(dataset.source)
        files_to_delete.update(dataset.pieces)

    deleted_count = 0
    for path in sorted(files_to_delete):
        if path.parent != input_dir:
            continue
        if path.suffix not in {".pvtu", ".vtu"}:
            continue
        if not path.exists():
            continue

        path.unlink()
        deleted_count += 1

    print(f"Deleted {deleted_count} parallel .pvtu/.vtu file(s) from {input_dir}.")


def merge_vtu_files(
    input_dir,
    expected_parts,
    simulation_name=None,
    overwrite=False,
    dry_run=False,
):
    input_dir = input_dir.resolve()
    datasets = collect_datasets(input_dir, expected_parts, simulation_name)
    if not datasets:
        print(f"No parallel VTU files found in {input_dir}.")
        return 1

    by_simulation = {}
    for dataset in datasets:
        by_simulation.setdefault(dataset.simulation, []).append(dataset)

    for simulation, simulation_datasets in sorted(by_simulation.items()):
        output_dir = input_dir / simulation
        if dry_run:
            print(
                f"Would process {simulation}: {len(simulation_datasets)} timestep(s) "
                f"into {output_dir}."
            )
            continue

        output_dir.mkdir(exist_ok=True)
        merged_datasets = []
        print(
            f"Processing {simulation}: {len(simulation_datasets)} timestep(s), "
            f"{simulation_datasets[0].num_parts} part(s)."
        )

        for dataset in sorted(simulation_datasets, key=lambda item: int(item.step)):
            output_path = output_dir / merged_file_name(dataset)

            if output_path.exists() and not overwrite:
                print(f"Skipping existing {output_path}")
                merged_datasets.append(dataset)
                continue

            if dataset.source is None and not dataset.pieces:
                print(f"Warning: no pieces found for {merged_file_name(dataset)}")
                continue

            merge_dataset(dataset, output_path)
            print(f"Merged {output_path}")
            merged_datasets.append(dataset)

        write_pvd(output_dir, simulation, merged_datasets)
        if len(merged_datasets) == len(simulation_datasets):
            delete_parallel_files(input_dir, merged_datasets)
        else:
            print("Not deleting parallel files because not all timesteps were merged.")

    return 0


def main():
    args = parse_args()
    if args.num_parts is not None and args.num_parts < 1:
        print("--num-parts must be greater than zero.", file=sys.stderr)
        return 2

    return merge_vtu_files(
        args.input_dir,
        args.num_parts,
        simulation_name=args.simulation_name,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    sys.exit(main())
