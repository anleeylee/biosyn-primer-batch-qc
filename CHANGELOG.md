# Changelog

## [1.0.0] - 2024-01-01
### Added
- Batch primer import (CSV / XLSX / TSV)
- Tm calculation: basic Wallace and nearest-neighbor (SantaLucia 1998)
- GC%, homopolymer, hairpin risk, self-complementarity checks
- Fully configurable QC rules (length, GC%, Tm ranges)
- Four-level QC: PASS / REVIEW / WARNING / INVALID
- Export to XLSX (multi-sheet), CSV, HTML
- Fault-tolerant batch processing (10,000+ primers)
- Offline / local-first operation
