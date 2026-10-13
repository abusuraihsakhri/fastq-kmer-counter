#!/usr/bin/env python3
"""
FASTQ K-mer Counter — Real Implementation

Parses FASTQ files, counts k-mers for arbitrary k, computes entropy,
quality score statistics (Phred), and sequence complexity estimation.

Zero external dependencies — Python stdlib only.
"""

import csv
import math
import os
import sys
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

# Maximum file size (100 MB) to prevent resource exhaustion
MAX_FILE_SIZE = 100 * 1024 * 1024


def _safe_resolve_path(filepath: str) -> Path:
    """Resolve and validate a file path, preventing path traversal."""
    p = Path(filepath).resolve()
    # Ensure path is within current working directory or an absolute allowed path
    cwd = Path.cwd().resolve()
    try:
        p.relative_to(cwd)
    except ValueError:
        # Allow paths under common temp directories for testing
        import tempfile
        temp_dir = Path(tempfile.gettempdir()).resolve()
        try:
            p.relative_to(temp_dir)
        except ValueError:
            raise ValueError(
                f"Path '{filepath}' is outside the allowed directory. "
                f"Paths must be within the project directory or system temp."
            )
    return p


# ---------------------------------------------------------------------------
# FASTQ Parsing
# ---------------------------------------------------------------------------

def parse_fastq(text: str) -> List[Dict[str, str]]:
    """Parse strict, four-line FASTQ records; reject malformed input.

    Accepts IUPAC nucleotide ambiguity codes. A quality character must be
    printable Phred-encoded ASCII and quality length must equal read length.
    """
    if not text.strip():
        return []
    lines = text.splitlines()
    if len(lines) % 4:
        raise ValueError("Incomplete FASTQ record: expected four lines per read.")

    records = []
    for i in range(0, len(lines), 4):
        header, sequence, plus, quality = (line.rstrip("\r") for line in lines[i:i + 4])
        record_number = i // 4 + 1
        if not header.startswith("@") or len(header) == 1:
            raise ValueError(f"FASTQ record {record_number}: expected an @header.")
        if not plus.startswith("+") or (len(plus) > 1 and plus[1:] != header[1:]):
            raise ValueError(f"FASTQ record {record_number}: invalid + separator.")
        if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence, re.IGNORECASE):
            raise ValueError(f"FASTQ record {record_number}: invalid nucleotide sequence.")
        if len(sequence) != len(quality):
            raise ValueError(f"FASTQ record {record_number}: quality length differs from sequence length.")
        if any(not 33 <= ord(c) <= 126 for c in quality):
            raise ValueError(f"FASTQ record {record_number}: invalid quality encoding.")
        records.append({
            "header": header[1:],
            "sequence": sequence.upper(),
            "quality": quality,
        })
    return records


def parse_fastq_file(filepath: str) -> List[Dict[str, str]]:
    """Parse a FASTQ file from disk."""
    safe_path = _safe_resolve_path(filepath)
    file_size = safe_path.stat().st_size
    if file_size > MAX_FILE_SIZE:
        raise ValueError(f"File too large ({file_size} bytes). Maximum allowed: {MAX_FILE_SIZE} bytes.")
    with open(safe_path, 'r', encoding='utf-8') as f:
        return parse_fastq(f.read())


# ---------------------------------------------------------------------------
# K-mer Counting
# ---------------------------------------------------------------------------

def count_kmers(sequences: List[str], k: int) -> Dict[str, int]:
    """Count all k-mers across a list of sequences.

    Args:
        sequences: List of DNA sequences.
        k: K-mer length.

    Returns:
        Dict mapping k-mer string to count.
    """
    if not isinstance(k, int) or k <= 0:
        return {}
    if k > 100:
        raise ValueError(f"k value {k} exceeds maximum allowed (100).")

    counts = Counter()
    for seq in sequences:
        seq = seq.upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if all(base in 'ACGT' for base in kmer):  # Exclude all IUPAC ambiguity codes
                counts[kmer] += 1

    return dict(counts)


def count_kmers_single(sequence: str, k: int) -> Dict[str, int]:
    """Count k-mers in a single sequence."""
    return count_kmers([sequence], k)


def kmer_frequency_spectrum(kmer_counts: Dict[str, int]) -> Dict[int, int]:
    """Compute frequency spectrum: histogram of k-mer counts.

    Returns dict mapping count_value -> number_of_kmers_with_that_count.
    """
    spectrum = Counter()
    for count in kmer_counts.values():
        spectrum[count] += 1
    return dict(sorted(spectrum.items()))


def most_frequent_kmers(kmer_counts: Dict[str, int],
                         n: int = 10) -> List[Tuple[str, int]]:
    """Return the N most frequent k-mers."""
    return sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)[:n]


def least_frequent_kmers(kmer_counts: Dict[str, int],
                          n: int = 10) -> List[Tuple[str, int]]:
    """Return the N least frequent k-mers."""
    return sorted(kmer_counts.items(), key=lambda x: x[1])[:n]


def unique_kmer_count(kmer_counts: Dict[str, int]) -> int:
    """Count number of distinct k-mers."""
    return len(kmer_counts)


def total_kmer_count(kmer_counts: Dict[str, int]) -> int:
    """Count total k-mer occurrences."""
    return sum(kmer_counts.values())


# ---------------------------------------------------------------------------
# Entropy Calculation
# ---------------------------------------------------------------------------

def kmer_entropy(kmer_counts: Dict[str, int]) -> float:
    """Calculate Shannon entropy of k-mer distribution.

    H = -Σ(p × log2(p))

    Higher entropy = more uniform distribution = more complex sequence.
    """
    total = sum(kmer_counts.values())
    if total == 0:
        return 0.0

    entropy = 0.0
    for count in kmer_counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)

    return entropy


def normalized_entropy(kmer_counts: Dict[str, int]) -> float:
    """Normalized entropy: H / log2(N) where N = number of distinct k-mers.

    Returns value between 0 and 1.
    """
    n = len(kmer_counts)
    if n <= 1:
        return 0.0

    h = kmer_entropy(kmer_counts)
    max_h = math.log2(n)
    return h / max_h if max_h > 0 else 0.0


# ---------------------------------------------------------------------------
# Quality Score Statistics (Phred Scores)
# ---------------------------------------------------------------------------

def phred_to_prob(phred: int) -> float:
    """Convert Phred score to error probability.

    P = 10^(-Q/10)
    """
    return 10.0 ** (-phred / 10.0)


def prob_to_phred(prob: float) -> float:
    """Convert error probability to Phred score.

    Q = -10 * log10(P)
    """
    if prob <= 0:
        return float('inf')
    return -10.0 * math.log10(prob)


def ascii_to_phred(char: str, offset: int = 33) -> int:
    """Convert ASCII character to Phred score.

    Standard Sanger/Illumina 1.8+ uses offset 33 (Phred+33).
    """
    return ord(char) - offset


def quality_scores(quality_string: str, offset: int = 33) -> List[int]:
    """Convert quality string to list of Phred scores."""
    return [ascii_to_phred(c, offset) for c in quality_string]


def per_read_quality(quality_string: str, offset: int = 33) -> Dict[str, Any]:
    """Calculate quality statistics for a single read."""
    scores = quality_scores(quality_string, offset)
    if not scores:
        return {'mean_quality': 0.0, 'min_quality': 0, 'max_quality': 0}

    mean_q = sum(scores) / len(scores)
    above_q20 = sum(1 for s in scores if s >= 20)
    above_q30 = sum(1 for s in scores if s >= 30)

    return {
        'mean_quality': round(mean_q, 2),
        'min_quality': min(scores),
        'max_quality': max(scores),
        'length': len(scores),
        'percent_above_q20': round(above_q20 / len(scores) * 100, 2),
        'percent_above_q30': round(above_q30 / len(scores) * 100, 2),
    }


def per_position_quality(records: List[Dict[str, str]],
                          offset: int = 33) -> List[Dict[str, Any]]:
    """Calculate mean quality at each position across all reads.

    Returns list of dicts with position, mean_quality, min_quality, max_quality.
    """
    if not records:
        return []

    max_len = max(len(r['quality']) for r in records)
    position_sums = [0.0] * max_len
    position_counts = [0] * max_len
    position_mins = [float('inf')] * max_len
    position_maxs = [0] * max_len

    for rec in records:
        scores = quality_scores(rec['quality'], offset)
        for i, s in enumerate(scores):
            position_sums[i] += s
            position_counts[i] += 1
            position_mins[i] = min(position_mins[i], s)
            position_maxs[i] = max(position_maxs[i], s)

    results = []
    for i in range(max_len):
        if position_counts[i] > 0:
            results.append({
                'position': i + 1,
                'mean_quality': round(position_sums[i] / position_counts[i], 2),
                'min_quality': position_mins[i],
                'max_quality': position_maxs[i],
                'count': position_counts[i],
            })

    return results


def overall_quality_stats(records: List[Dict[str, str]],
                           offset: int = 33) -> Dict[str, Any]:
    """Calculate overall quality statistics across all reads."""
    if not records:
        return {'num_reads': 0}

    all_scores = []
    per_read_means = []
    total_above_q20 = 0
    total_above_q30 = 0
    total_bases = 0

    for rec in records:
        scores = quality_scores(rec['quality'], offset)
        all_scores.extend(scores)
        per_read_means.append(sum(scores) / len(scores) if scores else 0)
        total_above_q20 += sum(1 for s in scores if s >= 20)
        total_above_q30 += sum(1 for s in scores if s >= 30)
        total_bases += len(scores)

    return {
        'num_reads': len(records),
        'total_bases': total_bases,
        'mean_quality': round(sum(all_scores) / len(all_scores), 2) if all_scores else 0,
        'min_quality': min(all_scores) if all_scores else 0,
        'max_quality': max(all_scores) if all_scores else 0,
        'mean_per_read_quality': round(sum(per_read_means) / len(per_read_means), 2),
        'percent_above_q20': round(total_above_q20 / total_bases * 100, 2) if total_bases else 0,
        'percent_above_q30': round(total_above_q30 / total_bases * 100, 2) if total_bases else 0,
        'read_length_min': min(len(r['sequence']) for r in records),
        'read_length_max': max(len(r['sequence']) for r in records),
    }


# ---------------------------------------------------------------------------
# Sequence Complexity
# ---------------------------------------------------------------------------

def sequence_complexity(sequence: str) -> Dict[str, Any]:
    """Estimate sequence complexity using multiple metrics.

    1. Shannon entropy of base composition
    2. Dinucleotide diversity
    3. Low-complexity run detection (homopolymer runs)
    """
    seq = sequence.upper()
    length = len(seq)
    if length == 0:
        return {'complexity': 0.0, 'entropy': 0.0}

    # 1. Base composition entropy
    base_counts = Counter(seq)
    base_entropy = 0.0
    for count in base_counts.values():
        if count > 0:
            p = count / length
            base_entropy -= p * math.log2(p)

    # Max entropy for 4 bases = log2(4) = 2.0
    normalized_base_entropy = base_entropy / 2.0

    # 2. Dinucleotide diversity
    dinuc_counts = Counter()
    for i in range(length - 1):
        dinuc_counts[seq[i:i + 2]] += 1
    n_dinucs = sum(dinuc_counts.values())
    dinuc_entropy = 0.0
    for count in dinuc_counts.values():
        if count > 0:
            p = count / n_dinucs
            dinuc_entropy -= p * math.log2(p)
    # Max entropy for 16 dinucleotides = log2(16) = 4.0
    normalized_dinuc_entropy = dinuc_entropy / 4.0 if n_dinucs > 0 else 0.0

    # 3. Homopolymer run detection
    max_run = 1
    current_run = 1
    for i in range(1, length):
        if seq[i] == seq[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1

    # Overall complexity: weighted average
    complexity = (normalized_base_entropy * 0.4 +
                  normalized_dinuc_entropy * 0.3 +
                  min(1.0, 1.0 - (max_run - 1) / length) * 0.3) if length > 1 else normalized_base_entropy

    return {
        'complexity': round(complexity, 4),
        'base_entropy': round(base_entropy, 4),
        'normalized_base_entropy': round(normalized_base_entropy, 4),
        'dinucleotide_entropy': round(dinuc_entropy, 4),
        'normalized_dinuc_entropy': round(normalized_dinuc_entropy, 4),
        'max_homopolymer_run': max_run,
        'length': length,
    }


# ---------------------------------------------------------------------------
# Full Analysis
# ---------------------------------------------------------------------------

def analyze_fastq(records: List[Dict[str, str]], k_values: List[int] = None,
                   offset: int = 33) -> Dict[str, Any]:
    """Complete analysis of FASTQ records.

    Args:
        records: Parsed FASTQ records.
        k_values: List of k values for k-mer analysis (default [1, 2, 3]).
        offset: Phred score offset (default 33).

    Returns:
        Comprehensive analysis dict.
    """
    if k_values is None:
        k_values = [1, 2, 3]
    # Validate k_values
    k_values = [k for k in k_values if isinstance(k, int) and 0 < k <= 100]
    if not k_values:
        k_values = [1, 2, 3]

    sequences = [r['sequence'] for r in records]

    # Quality stats
    quality = overall_quality_stats(records, offset)

    # K-mer analysis for each k
    kmer_results = {}
    for k in k_values:
        counts = count_kmers(sequences, k)
        kmer_results[f'k={k}'] = {
            'unique_kmers': unique_kmer_count(counts),
            'total_kmers': total_kmer_count(counts),
            'entropy': round(kmer_entropy(counts), 4),
            'normalized_entropy': round(normalized_entropy(counts), 4),
            'most_frequent': most_frequent_kmers(counts, 5),
            'least_frequent': least_frequent_kmers(counts, 5),
            'frequency_spectrum': kmer_frequency_spectrum(counts),
        }

    # Sequence complexity (average across reads)
    complexities = [sequence_complexity(r['sequence']) for r in records]
    avg_complexity = sum(c['complexity'] for c in complexities) / len(complexities) if complexities else 0

    return {
        'quality': quality,
        'kmer_analysis': kmer_results,
        'sequence_complexity': {
            'mean_complexity': round(avg_complexity, 4),
            'min_complexity': round(min(c['complexity'] for c in complexities), 4) if complexities else 0,
            'max_complexity': round(max(c['complexity'] for c in complexities), 4) if complexities else 0,
        },
    }


# ---------------------------------------------------------------------------
# Batch Processing
# ---------------------------------------------------------------------------

def process_batch(input_csv: str, output_csv: str,
                   k: int = 3) -> List[Dict[str, Any]]:
    """Process a CSV of sequences and count k-mers.

    Expected CSV columns: sequence (required).
    """
    if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= 100:
        raise ValueError("k must be an integer between 1 and 100.")
    safe_input = _safe_resolve_path(input_csv)
    safe_output = _safe_resolve_path(output_csv)
    if safe_input == safe_output:
        raise ValueError("Input and output CSV paths must differ.")
    if safe_input.stat().st_size > MAX_FILE_SIZE:
        raise ValueError("Input CSV exceeds the 100 MiB limit.")

    with open(safe_input, mode='r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        if not any(column in fieldnames for column in ('sequence', 'seq')):
            raise ValueError("CSV must contain a 'sequence' or 'seq' column.")
        rows = list(reader)

    out_fields = list(dict.fromkeys(fieldnames + [
        'length', 'gc_content', f'entropy_k{k}', 'complexity', 'num_unique_kmers'
    ]))
    out_rows = []

    for r in rows:
        seq = (r.get('sequence') or r.get('seq') or '').strip()
        if not seq:
            out_rows.append({**r, **{key: '' for key in ('length', 'gc_content', f'entropy_k{k}', 'complexity', 'num_unique_kmers')}})
            continue

        seq = seq.upper()
        counts = count_kmers_single(seq, k)
        ent = kmer_entropy(counts)
        comp = sequence_complexity(seq)

        # GC content
        gc = seq.count('G') + seq.count('C')
        atgc = seq.count('A') + seq.count('T') + gc
        gc_pct = gc / atgc * 100 if atgc > 0 else 0

        row_dict = dict(r)
        row_dict['length'] = len(seq)
        row_dict['gc_content'] = round(gc_pct, 2)
        row_dict[f'entropy_k{k}'] = round(ent, 4)
        row_dict['complexity'] = comp['complexity']
        row_dict['num_unique_kmers'] = len(counts)
        out_rows.append(row_dict)

    with open(safe_output, mode='w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    return out_rows
