# BioSyn Primer Batch QC

> **Local-first batch primer quality control** — check Tm, GC%, hairpin, and dimer risk for thousands of primers from a single CSV/XLSX file.

## Overview

Designing primers for a large experiment? Checking each one individually in a web tool is slow. **BioSyn Primer Batch QC** takes your primer list (CSV / XLSX / TSV) and runs a full computational QC sweep in seconds: Tm, GC%, homopolymer runs, hairpin risk, and self-complementarity — all with configurable thresholds.

> For single-primer online calculations, visit [BioSyn primer design tools](https://biosyn-inc.com/tools/primer-design-qc) — this is the batch-processing desktop companion.

## Quick Start

```bash
git clone https://github.com/your-org/biosyn-primer-batch-qc.git
cd biosyn-primer-batch-qc
pip install -r requirements.txt

python -m src.main examples/sample_primers.csv --output results.xlsx
python -m src.main primers.xlsx --tm-method nearest_neighbor --salt 50
```

## Features

- **Tm calculation**: Basic Wallace or Nearest-Neighbor (SantaLucia 1998)
- **GC% analysis** with configurable preferred range
- **Homopolymer detection** (AAAA, GGGG, CCCC, TTTT)
- **Hairpin risk**: Low / Review / High
- **Self-complementarity**: none / possible / strong
- **Fully configurable QC rules**
- **10,000+ primers** — single bad entry doesn't kill the run
- **CSV / XLSX / HTML export**
- **100% offline** — no primer sequences leave your machine

## Methods

| Metric | Method |
|--------|--------|
| Tm (basic) | 4×(G+C) + 2×(A+T), salt-corrected |
| Tm (NN) | SantaLucia 1998 unified nearest-neighbor parameters |
| Hairpin | Self-complementary stem-loop heuristic |
| Dimer | Overlapping reverse-complement stretch detection |

## Scientific Disclaimer

This is computational primer QC / decision support. It does **not** guarantee PCR success, absence of primer dimers, hairpins, or experimental specificity.

## Keywords

`primer batch QC`, `primer design tool`, `Tm calculator`, `hairpin prediction`, `primer dimer check`, `PCR primer QC`, `oligo analysis`, `batch primer validation`, `local bioinformatics`

## About

Built by [BioSyn Inc.](https://biosyn-inc.com/) — explore more [primer QC and PCR design tools](https://biosyn-inc.com/tools/primer-design-qc).

## License

[MIT](LICENSE)
