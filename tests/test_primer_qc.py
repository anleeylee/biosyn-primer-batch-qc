"""Tests for BioSyn Primer Batch QC."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from primer_analyzer import PrimerAnalyzer, QCConfig, evaluate_qc


def test_basic_primer():
    a = PrimerAnalyzer()
    r = a.analyze("P1", "ATCGATCGATCGATCGAT")
    assert r.length == 20
    assert r.gc_percent == 50.0
    assert r.valid


def test_gc_percent():
    a = PrimerAnalyzer()
    r = a.analyze("P2", "GCGCGCGCGC")
    assert r.gc_percent == 100.0


def test_homopolymer():
    a = PrimerAnalyzer()
    r = a.analyze("P3", "ATGGGGATCG")
    assert "GGGG" in r.homopolymer


def test_qc_pass():
    a = PrimerAnalyzer(QCConfig())
    r = a.analyze("P4", "ATCGATCGATCGATCGAT")
    qc = evaluate_qc(r, QCConfig())
    assert qc.status in ("PASS", "REVIEW")


def test_invalid_sequence():
    a = PrimerAnalyzer()
    r = a.analyze("P5", "ATCGXYZ")
    assert not r.valid


def test_empty():
    a = PrimerAnalyzer()
    r = a.analyze("P6", "")
    assert not r.valid


if __name__ == "__main__":
    test_basic_primer(); test_gc_percent(); test_homopolymer()
    test_qc_pass(); test_invalid_sequence(); test_empty()
    print("All tests passed!")
