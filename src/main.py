"""
CLI entry point for BioSyn Primer Batch QC.

Usage:
    python -m src.main input.xlsx --output results.xlsx
    python -m src.main primers.csv --tm-method nearest_neighbor --salt 50
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any

from .primer_analyzer import (
    PrimerAnalyzer, QCConfig, PrimerResult,
    evaluate_qc, PrimerQCStatus,
)
from .exporter import export_primer_results

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def read_primer_file(path: Path) -> List[Dict[str, str]]:
    """
    Read primer CSV / XLSX / TSV file.

    Required columns: Primer_ID, Sequence
    Optional: Name, Target, Gene, Direction, Sample, Notes
    """
    records: List[Dict[str, str]] = []

    if path.suffix.lower() in (".xlsx", ".xls"):
        try:
            import openpyxl
            wb = openpyxl.load_workbook(path, read_only=True)
            ws = wb.active
            rows = list(ws.iter_rows(values_only=True))
            if not rows:
                return records
            headers = [str(h).strip() if h else "" for h in rows[0]]
            for row in rows[1:]:
                record = {}
                for h, v in zip(headers, row):
                    record[h] = str(v).strip() if v else ""
                records.append(record)
        except ImportError:
            raise RuntimeError("openpyxl required for XLSX files. pip install openpyxl")
    else:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter=delimiter)
            for row in reader:
                records.append({k: (v or "").strip() for k, v in row.items()})

    logger.info("Read %d primer records from %s", len(records), path)
    return records


def run_batch_qc(
    input_path: Path,
    output_path: Path,
    config: QCConfig,
) -> Dict[str, Any]:
    """Run the full batch primer QC pipeline."""
    # Read input
    records = read_primer_file(input_path)
    if not records:
        raise ValueError("No primer records found in input file.")

    # Map column names (flexible)
    def find_key(keys: List[str]) -> str:
        """Find a matching column key case-insensitively."""
        for k in keys:
            for actual in records[0].keys():
                if actual.lower() == k.lower():
                    return actual
        return keys[0]  # fallback

    id_key = find_key(["Primer_ID", "ID", "primer_id", "id", "Name"])
    seq_key = find_key(["Sequence", "Seq", "sequence", "seq"])

    # Analyze each primer independently
    analyzer = PrimerAnalyzer(config)
    results: List[PrimerResult] = []
    qc_statuses: List[PrimerQCStatus] = []
    errors: List[Dict[str, str]] = []

    for i, rec in enumerate(records):
        pid = rec.get(id_key, f"primer_{i+1}")
        seq = rec.get(seq_key, "")
        try:
            result = analyzer.analyze(pid, seq)
            results.append(result)
            qc_statuses.append(evaluate_qc(result, config))
        except Exception as exc:
            logger.error("Failed to analyze %s: %s", pid, exc)
            results.append(PrimerResult(
                primer_id=pid, sequence=seq, valid=False, error=str(exc),
            ))
            qc_statuses.append(PrimerQCStatus(
                primer_id=pid, status="INVALID", warnings=[str(exc)],
            ))
            errors.append({"primer_id": pid, "error": str(exc)})

    # Export
    output_path.parent.mkdir(parents=True, exist_ok=True)
    export_primer_results(results, qc_statuses, config, errors, output_path)

    summary = {
        "total": len(results),
        "PASS": sum(1 for q in qc_statuses if q.status == "PASS"),
        "REVIEW": sum(1 for q in qc_statuses if q.status == "REVIEW"),
        "WARNING": sum(1 for q in qc_statuses if q.status == "WARNING"),
        "INVALID": sum(1 for q in qc_statuses if q.status == "INVALID"),
    }
    logger.info("Batch QC complete: %s", summary)
    return summary


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="biosyn-primer-qc",
        description="BioSyn Primer Batch QC — batch primer quality control.",
        epilog="Built by BioSyn Inc. — https://biosyn-inc.com/",
    )
    parser.add_argument("input", type=Path, help="Input CSV / XLSX / TSV file")
    parser.add_argument("--output", "-o", type=Path,
                        default=Path("primer_qc_results.xlsx"),
                        help="Output file path")
    parser.add_argument("--tm-method", choices=["basic", "nearest_neighbor"],
                        default="basic", help="Tm calculation method")
    parser.add_argument("--salt", type=float, default=50.0,
                        help="Salt concentration (mM)")
    parser.add_argument("--mg", type=float, default=1.5,
                        help="Mg2+ concentration (mM)")
    parser.add_argument("--primer-conc", type=float, default=0.25,
                        help="Primer concentration (µM)")
    parser.add_argument("--min-length", type=int, default=18)
    parser.add_argument("--max-length", type=int, default=25)
    parser.add_argument("--gc-min", type=float, default=40.0)
    parser.add_argument("--gc-max", type=float, default=60.0)
    parser.add_argument("--tm-min", type=float, default=55.0)
    parser.add_argument("--tm-max", type=float, default=65.0)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    setup_logging(args.verbose)

    config = QCConfig(
        min_length=args.min_length,
        max_length=args.max_length,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        tm_min=args.tm_min,
        tm_max=args.tm_max,
        salt_conc=args.salt,
        mg_conc=args.mg,
        primer_conc=args.primer_conc,
        tm_method=args.tm_method,
    )

    try:
        summary = run_batch_qc(args.input, args.output, config)
        print(f"\nBatch QC complete: {summary}")
        print(f"Results saved to: {args.output.resolve()}")
        return 0
    except Exception as exc:
        logger.error("Fatal error: %s", exc)
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
