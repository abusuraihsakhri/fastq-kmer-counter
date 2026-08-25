#!/usr/bin/env python3
"""Tests for FASTQ K-mer Counter — 20+ real functional tests."""
import pytest
from kmer_counter import (
    parse_fastq,
    count_kmers,
    count_kmers_single,
    kmer_frequency_spectrum,
    most_frequent_kmers,
    least_frequent_kmers,
    unique_kmer_count,
    total_kmer_count,
    kmer_entropy,
    normalized_entropy,
    phred_to_prob,
    prob_to_phred,
    ascii_to_phred,
    quality_scores,
    per_read_quality,
    overall_quality_stats,
    sequence_complexity,
    analyze_fastq,
    process_batch,
)


# ---------------------------------------------------------------------------
# FASTQ Parsing Tests
# ---------------------------------------------------------------------------

class TestFASTQParsing:
    def test_parse_single_record(self):
        text = "@read1\nATCGATCG\n+\nIIIIIIII\n"
        records = parse_fastq(text)
        assert len(records) == 1
        assert records[0]['header'] == 'read1'
        assert records[0]['sequence'] == 'ATCGATCG'
        assert records[0]['quality'] == 'IIIIIIII'

    def test_parse_multiple_records(self):
        text = "@r1\nATCG\n+\nIIII\n@r2\nGGCC\n+\nIIII\n"
        records = parse_fastq(text)
        assert len(records) == 2

    def test_parse_empty(self):
        records = parse_fastq('')
        assert len(records) == 0

    def test_sequence_uppercased(self):
        text = "@r1\natcg\n+\nIIII\n"
        records = parse_fastq(text)
        assert records[0]['sequence'] == 'ATCG'


# ---------------------------------------------------------------------------
# K-mer Counting Tests
# ---------------------------------------------------------------------------

class TestKmerCounting:
    def test_1mer(self):
        counts = count_kmers_single('ATCG', 1)
        assert counts['A'] == 1
        assert counts['T'] == 1
        assert counts['C'] == 1
        assert counts['G'] == 1

    def test_2mer(self):
        counts = count_kmers_single('ATCG', 2)
        assert counts['AT'] == 1
        assert counts['TC'] == 1
        assert counts['CG'] == 1

    def test_3mer(self):
        counts = count_kmers_single('ATCGATCG', 3)
        assert counts['ATC'] == 2
        assert counts['TCG'] == 2
        assert counts['CGA'] == 1

    def test_repeated_sequence(self):
        counts = count_kmers_single('AAAA', 2)
        assert counts['AA'] == 3

    def test_skip_n_bases(self):
        counts = count_kmers_single('ANCG', 2)
        assert 'AN' not in counts
        assert 'NC' not in counts
        assert 'CG' in counts

    def test_k_too_large(self):
        counts = count_kmers_single('AT', 5)
        assert len(counts) == 0

    def test_multiple_sequences(self):
        counts = count_kmers(['ATCG', 'ATCG'], 2)
        assert counts['AT'] == 2


# ---------------------------------------------------------------------------
# Frequency Spectrum Tests
# ---------------------------------------------------------------------------

class TestFrequencySpectrum:
    def test_basic_spectrum(self):
        counts = {'AA': 3, 'AT': 3, 'CG': 1}
        spectrum = kmer_frequency_spectrum(counts)
        assert spectrum[1] == 1  # One k-mer with count 1
        assert spectrum[3] == 2  # Two k-mers with count 3

    def test_empty(self):
        spectrum = kmer_frequency_spectrum({})
        assert len(spectrum) == 0


# ---------------------------------------------------------------------------
# Most/Least Frequent Tests
# ---------------------------------------------------------------------------

class TestFrequentKmers:
    def test_most_frequent(self):
        counts = {'AA': 5, 'AT': 3, 'CG': 1}
        top = most_frequent_kmers(counts, 2)
        assert top[0] == ('AA', 5)
        assert top[1] == ('AT', 3)

    def test_least_frequent(self):
        counts = {'AA': 5, 'AT': 3, 'CG': 1}
        bottom = least_frequent_kmers(counts, 1)
        assert bottom[0] == ('CG', 1)


# ---------------------------------------------------------------------------
# Entropy Tests
# ---------------------------------------------------------------------------

class TestEntropy:
    def test_uniform_distribution(self):
        # All k-mers equally frequent -> maximum entropy
        counts = {'A': 10, 'T': 10, 'G': 10, 'C': 10}
        ent = kmer_entropy(counts)
        assert abs(ent - 2.0) < 0.01  # log2(4) = 2.0

    def test_single_kmer(self):
        # Only one k-mer -> zero entropy
        counts = {'A': 100}
        ent = kmer_entropy(counts)
        assert ent == 0.0

    def test_empty(self):
        assert kmer_entropy({}) == 0.0

    def test_normalized(self):
        counts = {'A': 10, 'T': 10, 'G': 10, 'C': 10}
        norm = normalized_entropy(counts)
        assert abs(norm - 1.0) < 0.01  # Perfectly uniform

    def test_normalized_single(self):
        counts = {'A': 100}
        assert normalized_entropy(counts) == 0.0


# ---------------------------------------------------------------------------
# Quality Score Tests
# ---------------------------------------------------------------------------

class TestQualityScores:
    def test_phred_to_prob(self):
        # Q30 = 0.001 error rate
        assert abs(phred_to_prob(30) - 0.001) < 1e-6

    def test_phred_to_prob_q20(self):
        assert abs(phred_to_prob(20) - 0.01) < 1e-6

    def test_prob_to_phred(self):
        assert abs(prob_to_phred(0.001) - 30.0) < 0.01

    def test_ascii_to_phred(self):
        # 'I' = ASCII 73, Phred+33 -> 73-33 = 40
        assert ascii_to_phred('I') == 40

    def test_quality_scores(self):
        scores = quality_scores('IIII')
        assert all(s == 40 for s in scores)

    def test_per_read_quality(self):
        result = per_read_quality('IIIIIIII')
        assert result['mean_quality'] == 40.0
        assert result['percent_above_q30'] == 100.0

    def test_per_read_quality_low(self):
        # '!' = ASCII 33 -> Phred 0
        result = per_read_quality('!!!!')
        assert result['mean_quality'] == 0.0
        assert result['percent_above_q20'] == 0.0


# ---------------------------------------------------------------------------
# Overall Quality Stats Tests
# ---------------------------------------------------------------------------

class TestOverallQuality:
    def test_basic_stats(self):
        records = [
            {'sequence': 'ATCG', 'quality': 'IIII'},
            {'sequence': 'GGCC', 'quality': '!!!!'},
        ]
        stats = overall_quality_stats(records)
        assert stats['num_reads'] == 2
        assert stats['total_bases'] == 8

    def test_empty(self):
        stats = overall_quality_stats([])
        assert stats['num_reads'] == 0


# ---------------------------------------------------------------------------
# Sequence Complexity Tests
# ---------------------------------------------------------------------------

class TestSequenceComplexity:
    def test_high_complexity(self):
        result = sequence_complexity('ATCGATCGATCG')
        assert result['complexity'] > 0.8

    def test_low_complexity(self):
        result = sequence_complexity('AAAAAAAAAA')
        assert result['complexity'] < 0.3

    def test_empty(self):
        result = sequence_complexity('')
        assert result['complexity'] == 0.0

    def test_max_homopolymer(self):
        result = sequence_complexity('AAATCG')
        assert result['max_homopolymer_run'] == 3

    def test_components_present(self):
        result = sequence_complexity('ATCGATCG')
        assert 'base_entropy' in result
        assert 'dinucleotide_entropy' in result
        assert 'max_homopolymer_run' in result


# ---------------------------------------------------------------------------
# Full Analysis Tests
# ---------------------------------------------------------------------------

class TestFullAnalysis:
    def test_analyze(self):
        records = [
            {'header': 'r1', 'sequence': 'ATCGATCG', 'quality': 'IIIIIIII'},
        ]
        result = analyze_fastq(records, [1, 2, 3])
        assert 'quality' in result
        assert 'kmer_analysis' in result
        assert 'k=1' in result['kmer_analysis']
        assert 'k=3' in result['kmer_analysis']


# ---------------------------------------------------------------------------
# Batch Processing Tests
# ---------------------------------------------------------------------------

class TestBatchProcessing:
    def test_batch_csv(self, tmp_path):
        csv_in = tmp_path / 'in.csv'
        csv_out = tmp_path / 'out.csv'
        csv_in.write_text(
            'sequence\n'
            'ATCGATCGATCG\n'
            'AAAAAAAAAA\n',
            encoding='utf-8'
        )
        results = process_batch(str(csv_in), str(csv_out), k=3)
        assert len(results) == 2
        assert csv_out.exists()
        assert 'entropy_k3' in results[0]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
