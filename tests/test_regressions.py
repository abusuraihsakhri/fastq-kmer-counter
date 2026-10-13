"""Regression tests for real FASTQ analytics, CSV processing and audit security."""
import csv

import pytest

from kmer_counter import parse_fastq, count_kmers_single, process_batch
from agents.base import AuditTrail


@pytest.mark.parametrize("text", [
    "@read1\nACGT\n+\nIII\n",
    "@read1\nACGT\n+\nIIII\n@orphan\n",
    "@read1\nACGT\n-\nIIII\n",
    "@read1\nAC1T\n+\nIIII\n",
    "@read1\nACGT\n+other\nIIII\n",
    "@read1\nACGT\n+\nIIII\n\n",
])
def test_malformed_fastq_is_not_silently_accepted(text):
    with pytest.raises(ValueError):
        parse_fastq(text)


def test_ambiguous_kmers_are_excluded():
    assert count_kmers_single("ACRGT", 2) == {"AC": 1, "GT": 1}
    assert count_kmers_single("ACN", 1) == {"A": 1, "C": 1}


def test_batch_accepts_seq_and_deduplicates_output_fields(tmp_path):
    source = tmp_path / "in.csv"
    output = tmp_path / "out.csv"
    source.write_text("seq,entropy_k2\nACGT,old\n,unused\n", encoding="utf-8")
    result = process_batch(str(source), str(output), k=2)
    assert len(result) == 2
    assert result[0]["num_unique_kmers"] == 3
    with output.open(newline="", encoding="utf-8") as f:
        rows = csv.DictReader(f)
        assert rows.fieldnames.count("entropy_k2") == 1
        assert len(list(rows)) == 2


def test_batch_validates_header_k_and_paths(tmp_path):
    source = tmp_path / "in.csv"
    output = tmp_path / "out.csv"
    source.write_text("other\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="sequence.*seq"):
        process_batch(str(source), str(output), k=2)
    source.write_text("sequence\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="k must"):
        process_batch(str(source), str(output), k=0)
    with pytest.raises(ValueError, match="must differ"):
        process_batch(str(source), str(source), k=2)


def test_audit_recalculates_hmac_and_does_not_expose_mutable_logs():
    audit = AuditTrail(secret_key="test-only-secret")
    audit.log("worker", "analysis", "event", {"status": "ok"})
    assert audit.verify_integrity()
    audit.get_trail()[0]["current_hash"] = "injected"
    assert audit.verify_integrity()
    audit.logs[0]["actor"] = "tampered"
    assert not audit.verify_integrity()


def test_audit_rejects_consistently_relinked_forged_hashes():
    audit = AuditTrail(secret_key="test-only-secret")
    audit.log("worker", "analysis", "event", {"status": "ok"})
    audit.log("worker", "analysis", "event2", {"status": "ok"})
    audit.logs[0]["current_hash"] = "0" * 64
    audit.logs[1]["prev_hash"] = "0" * 64
    assert not audit.verify_integrity()
