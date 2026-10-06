"""
Primer analysis engine — Tm, GC%, homopolymer, hairpin, self-complementarity.

All calculations are computational QC estimates, not experimental guarantees.
"""

from __future__ import annotations

import re
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# DNA bases
VALID_BASES = set("ACGT")
AMBIG_BASES = set("RYSWKMBDHVN")

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998, unified)
# ΔH (kcal/mol), ΔS (cal/mol·K)
NN_PARAMS: Dict[str, Tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}


@dataclass
class PrimerResult:
    """Analysis result for a single primer."""

    primer_id: str
    sequence: str
    length: int = 0
    gc_percent: float = 0.0
    tm: float = 0.0
    tm_method: str = ""
    base_counts: Dict[str, int] = field(default_factory=dict)
    ambiguous_count: int = 0
    homopolymer: str = ""          # e.g. "AAAA(4)"
    repeat_warning: bool = False
    gc_clamp: bool = False          # 3' end ends with G or C
    three_prime_gc: int = 0         # number of G/C in last 3 bases
    self_complementarity: str = "none"  # none | possible | strong
    hairpin_risk: str = "Low"       # Low | Review | High
    valid: bool = True
    error: str = ""


@dataclass
class QCConfig:
    """Configurable QC thresholds."""

    min_length: int = 18
    max_length: int = 25
    gc_min: float = 40.0
    gc_max: float = 60.0
    tm_min: float = 55.0
    tm_max: float = 65.0
    homopolymer_max: int = 4        # max consecutive same base
    salt_conc: float = 50.0         # mM
    mg_conc: float = 1.5            # mM
    primer_conc: float = 0.25      # µM
    tm_method: str = "basic"        # basic | nearest_neighbor


class PrimerAnalyzer:
    """Core primer QC analysis engine."""

    def __init__(self, config: Optional[QCConfig] = None) -> None:
        self.config = config or QCConfig()

    def analyze(self, primer_id: str, sequence: str) -> PrimerResult:
        """
        Run full QC analysis on a single primer sequence.

        Args:
            primer_id: Unique primer identifier.
            sequence: Primer sequence (DNA bases).

        Returns:
            PrimerResult with all computed metrics.
        """
        seq = sequence.strip().upper() if sequence else ""
        result = PrimerResult(primer_id=primer_id, sequence=seq)

        if not seq:
            result.valid = False
            result.error = "Empty sequence"
            return result

        result.length = len(seq)
        result.base_counts = {b: seq.count(b) for b in "ACGTN"}
        result.ambiguous_count = sum(
            seq.count(b) for b in (AMBIG_BASES - {"N"})
        ) + seq.count("N")

        # Check validity
        invalid = set(seq) - (VALID_BASES | AMBIG_BASES)
        if invalid:
            result.valid = False
            result.error = f"Invalid characters: {''.join(sorted(invalid))}"

        # GC%
        total = result.base_counts["A"] + result.base_counts["C"] + \
                result.base_counts["G"] + result.base_counts["T"]
        if total > 0:
            result.gc_percent = round(
                (result.base_counts["G"] + result.base_counts["C"]) / total * 100, 2
            )

        # Tm
        result.tm_method = self.config.tm_method
        result.tm = round(self._compute_tm(seq), 2)

        # Homopolymer
        result.homopolymer = self._check_homopolymer(seq)

        # GC clamp (last base is G/C)
        result.gc_clamp = seq[-1] in "GC"
        result.three_prime_gc = sum(1 for b in seq[-3:] if b in "GC")

        # Self-complementarity
        result.self_complementarity = self._check_self_complementarity(seq)

        # Hairpin risk
        result.hairpin_risk = self._check_hairpin(seq)

        return result

    def _compute_tm(self, seq: str) -> float:
        """Compute melting temperature using configured method."""
        if self.config.tm_method == "nearest_neighbor":
            return self._tm_nearest_neighbor(seq)
        return self._tm_basic(seq)

    def _tm_basic(self, seq: str) -> float:
        """
        Basic approximation: Wallace / Marmur formula.
        Tm = 4*(G+C) + 2*(A+T) — adjusted for salt.
        """
        gc = seq.count("G") + seq.count("C")
        at = seq.count("A") + seq.count("T")
        if gc + at == 0:
            return 0.0
        tm = 4.0 * gc + 2.0 * at
        # Salt correction (rough approximation)
        na = self.config.salt_conc / 1000.0  # convert mM to M
        if na > 0:
            tm += 16.6 * math.log10(na)
        return tm

    def _tm_nearest_neighbor(self, seq: str) -> float:
        """
        Nearest-neighbor Tm (SantaLucia 1998 unified parameters).
        Tm = ΔH / (ΔS + R * ln(CT/4)) — 273.15 → °C
        """
        valid_seq = "".join(b for b in seq if b in "ACGT")
        if len(valid_seq) < 2:
            return self._tm_basic(valid_seq)

        dh_total = 0.0
        ds_total = 0.0
        for i in range(len(valid_seq) - 1):
            pair = valid_seq[i:i+2]
            if pair in NN_PARAMS:
                dh, ds = NN_PARAMS[pair]
                dh_total += dh
                ds_total += ds

        # Initiation penalties
        dh_total += 0.2   # initiation w/ terminal AT
        ds_total += -5.7

        r = 1.9872  # cal/mol·K
        ct = self.config.primer_conc * 1e-6  # convert µM to M
        # Salt correction for Mg
        mg = self.config.mg_conc * 1e-3
        salt_factor = 1.0 + 0.36 * math.sqrt(self.config.salt_conc / 1000.0)
        ds_total = ds_total * salt_factor + 0.36 * len(valid_seq) * math.log(
            self.config.salt_conc / 1000.0
        )

        try:
            tm_k = (dh_total * 1000) / (ds_total + r * math.log(ct / 4))
            return tm_k - 273.15
        except (ZeroDivisionError, ValueError):
            return self._tm_basic(valid_seq)

    def _check_homopolymer(self, seq: str) -> str:
        """Detect longest homopolymer run."""
        max_run = 1
        max_base = seq[0]
        current_run = 1
        for i in range(1, len(seq)):
            if seq[i] == seq[i-1]:
                current_run += 1
                if current_run > max_run:
                    max_run = current_run
                    max_base = seq[i]
            else:
                current_run = 1
        if max_run >= self.config.homopolymer_max:
            return f"{max_base * max_run}({max_run})"
        return ""

    def _check_self_complementarity(self, seq: str) -> str:
        """
        Check for self-complementary regions (primer-dimer risk indicator).
        Returns: none | possible | strong
        """
        comp = str.maketrans("ACGT", "TGCA")
        rev_comp = seq.translate(comp)[::-1]

        # Check for overlapping complementary stretches
        best_match = 0
        for k in range(4, min(len(seq) // 2, 10) + 1):
            for i in range(len(seq) - k + 1):
                subseq = seq[i:i+k]
                if subseq in rev_comp:
                    best_match = max(best_match, k)

        if best_match >= 6:
            return "strong"
        if best_match >= 4:
            return "possible"
        return "none"

    def _check_hairpin(self, seq: str) -> str:
        """
        Rough hairpin risk assessment based on self-complementary stems.
        Returns: Low | Review | High
        """
        comp = str.maketrans("ACGT", "TGCA")
        score = 0

        # Look for palindromic-ish regions
        for k in range(3, 8):
            for i in range(len(seq) - 2*k + 1):
                left = seq[i:i+k]
                right = seq[i+k:i+2*k]
                right_rc = right.translate(comp)[::-1]
                matches = sum(1 for a, b in zip(left, right_rc) if a == b)
                if matches >= k - 1:
                    score += 1

        if score >= 3:
            return "High"
        if score >= 1:
            return "Review"
        return "Low"


@dataclass
class PrimerQCStatus:
    """Final QC status for a primer."""

    primer_id: str
    status: str  # PASS | REVIEW | WARNING | INVALID
    warnings: List[str] = field(default_factory=list)


def evaluate_qc(result: PrimerResult, config: QCConfig) -> PrimerQCStatus:
    """
    Evaluate QC rules on a primer analysis result.

    Args:
        result: Primer analysis result.
        config: QC configuration with thresholds.

    Returns:
        PrimerQCStatus with PASS / REVIEW / WARNING / INVALID.
    """
    warnings: List[str] = []
    status = "PASS"

    if not result.valid:
        return PrimerQCStatus(
            primer_id=result.primer_id, status="INVALID",
            warnings=[result.error],
        )

    # Length
    if result.length < config.min_length or result.length > config.max_length:
        warnings.append(
            f"Length {result.length} outside preferred range "
            f"[{config.min_length}-{config.max_length}]"
        )
        if status == "PASS":
            status = "REVIEW"

    # GC%
    if result.gc_percent < config.gc_min or result.gc_percent > config.gc_max:
        warnings.append(
            f"GC% {result.gc_percent:.1f} outside preferred range "
            f"[{config.gc_min}-{config.gc_max}]"
        )
        if status == "PASS":
            status = "REVIEW"

    # Tm
    if result.tm < config.tm_min or result.tm > config.tm_max:
        warnings.append(
            f"Tm {result.tm:.1f}°C outside preferred range "
            f"[{config.tm_min}-{config.tm_max}]°C"
        )
        if status == "PASS":
            status = "REVIEW"

    # Homopolymer
    if result.homopolymer:
        warnings.append(f"Homopolymer: {result.homopolymer}")
        status = "WARNING"

    # Self-complementarity
    if result.self_complementarity == "strong":
        warnings.append("Strong self-complementarity — potential dimer candidate")
        status = "WARNING"
    elif result.self_complementarity == "possible":
        warnings.append("Possible self-complementarity — review")
        if status == "PASS":
            status = "REVIEW"

    # Hairpin
    if result.hairpin_risk == "High":
        warnings.append("High hairpin risk detected")
        status = "WARNING"
    elif result.hairpin_risk == "Review":
        warnings.append("Possible hairpin — review")
        if status == "PASS":
            status = "REVIEW"

    # Ambiguity
    if result.ambiguous_count > 0:
        warnings.append(f"{result.ambiguous_count} ambiguous base(s)")
        if status == "PASS":
            status = "REVIEW"

    return PrimerQCStatus(
        primer_id=result.primer_id, status=status, warnings=warnings,
    )
