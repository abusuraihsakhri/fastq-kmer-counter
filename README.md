# FASTQ K-mer Counter

Count k-mers, compute entropy, quality statistics, and sequence complexity from FASTQ files.

## Features

- **K-mer Counting**: Arbitrary k values, most/least frequent k-mers
- **Frequency Spectrum**: Histogram of k-mer count distributions
- **Shannon Entropy**: H = -Σ(p × log2(p)), normalized entropy
- **Quality Statistics**: Phred scores, per-position and per-read quality
- **Q20/Q30**: Percentage of bases above quality thresholds
- **Sequence Complexity**: Base entropy, dinucleotide diversity, homopolymer detection

## Installation

Zero dependencies — Python 3.8+ stdlib only.

```bash
cd fastq-kmer-counter
```

## Usage

### CLI

```bash
# Count 3-mers from a FASTQ file
python cli.py count -f reads.fastq -k 3

# Count k-mers from a sequence
python cli.py count -s "ATCGATCGATCG" -k 3

# Quality statistics
python cli.py quality -f reads.fastq

# Shannon entropy
python cli.py entropy -f reads.fastq -k 5

# Sequence complexity
python cli.py complexity -s "ATCGATCGATCG"

# Full analysis
python cli.py analyze -f reads.fastq --k-values 1,2,3,5

# Batch process from CSV
python cli.py batch -i sequences.csv -o results.csv -k 3
```

### Python API

```python
from kmer_counter import (
    count_kmers_single,
    kmer_entropy,
    sequence_complexity,
    per_read_quality,
    parse_fastq,
)

# Count k-mers
counts = count_kmers_single('ATCGATCGATCG', 3)

# Entropy
h = kmer_entropy(counts)

# Complexity
comp = sequence_complexity('ATCGATCGATCG')
print(f"Complexity: {comp['complexity']}")

# Quality
qual = per_read_quality('IIIIIIIIII')
print(f"Mean Q: {qual['mean_quality']}, >Q30: {qual['percent_above_q30']}%")
```

## Formulas

### Shannon Entropy
H = -Σ(p × log2(p)) where p = kmer_count / total_kmers

### Phred Score
Q = -10 × log10(P_error)
- Q20 = 1% error rate
- Q30 = 0.1% error rate

### Sequence Complexity
Weighted combination of:
- Base composition entropy (40%)
- Dinucleotide entropy (30%)
- Homopolymer run penalty (30%)

## Tests

```bash
python -m pytest test_kmer_counter.py -v
```

## License

MIT
