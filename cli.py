#!/usr/bin/env python3
"""
CLI for FASTQ K-mer Counter.

Subcommands:
  count     — Count k-mers in sequences
  quality   — Quality score statistics
  entropy   — Shannon entropy of k-mer distribution
  complexity — Sequence complexity estimation
  analyze   — Full FASTQ analysis
  batch     — Batch process from CSV
"""
import argparse
import json
import sys

from kmer_counter import (
    count_kmers,
    count_kmers_single,
    kmer_frequency_spectrum,
    most_frequent_kmers,
    least_frequent_kmers,
    kmer_entropy,
    normalized_entropy,
    per_read_quality,
    per_position_quality,
    overall_quality_stats,
    sequence_complexity,
    analyze_fastq,
    parse_fastq,
    parse_fastq_file,
    process_batch,
)


def cmd_count(args):
    """Count k-mers."""
    if args.file:
        records = parse_fastq_file(args.file)
        sequences = [r['sequence'] for r in records]
    elif args.sequence:
        sequences = [args.sequence]
    else:
        print("Provide --file or --sequence", file=sys.stderr)
        sys.exit(1)

    counts = count_kmers(sequences, args.k)
    top = most_frequent_kmers(counts, args.top)

    result = {
        'k': args.k,
        'unique_kmers': len(counts),
        'total_kmers': sum(counts.values()),
        'top_kmers': top,
    }
    print(json.dumps(result, indent=2))


def cmd_quality(args):
    """Quality score statistics."""
    if args.file:
        records = parse_fastq_file(args.file)
    else:
        print("Provide --file", file=sys.stderr)
        sys.exit(1)

    stats = overall_quality_stats(records)
    print(json.dumps(stats, indent=2))


def cmd_entropy(args):
    """Shannon entropy of k-mer distribution."""
    if args.file:
        records = parse_fastq_file(args.file)
        sequences = [r['sequence'] for r in records]
    elif args.sequence:
        sequences = [args.sequence]
    else:
        print("Provide --file or --sequence", file=sys.stderr)
        sys.exit(1)

    counts = count_kmers(sequences, args.k)
    result = {
        'k': args.k,
        'entropy': round(kmer_entropy(counts), 4),
        'normalized_entropy': round(normalized_entropy(counts), 4),
        'unique_kmers': len(counts),
    }
    print(json.dumps(result, indent=2))


def cmd_complexity(args):
    """Sequence complexity estimation."""
    if args.sequence:
        result = sequence_complexity(args.sequence)
        print(json.dumps(result, indent=2))
    elif args.file:
        records = parse_fastq_file(args.file)
        results = []
        for r in records[:10]:  # First 10 reads
            comp = sequence_complexity(r['sequence'])
            comp['header'] = r['header']
            results.append(comp)
        print(json.dumps(results, indent=2))
    else:
        print("Provide --sequence or --file", file=sys.stderr)
        sys.exit(1)


def cmd_analyze(args):
    """Full FASTQ analysis."""
    if args.file:
        records = parse_fastq_file(args.file)
    else:
        print("Provide --file", file=sys.stderr)
        sys.exit(1)

    k_values = [int(x) for x in args.k_values.split(',')]
    result = analyze_fastq(records, k_values)
    print(json.dumps(result, indent=2))


def cmd_batch(args):
    """Batch process from CSV."""
    results = process_batch(args.input, args.output, args.k)
    print(f"Processed {len(results)} sequences -> {args.output}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='fastq-kmer-counter',
        description='FASTQ K-mer Counter — K-mer counting, entropy, quality stats, complexity'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # count
    p_count = subparsers.add_parser('count', help='Count k-mers')
    p_count.add_argument('--sequence', '-s', help='DNA sequence')
    p_count.add_argument('--file', '-f', help='FASTQ file')
    p_count.add_argument('-k', type=int, default=3, help='K-mer size (default: 3)')
    p_count.add_argument('--top', type=int, default=10, help='Show top N k-mers')
    p_count.set_defaults(func=cmd_count)

    # quality
    p_qual = subparsers.add_parser('quality', help='Quality score statistics')
    p_qual.add_argument('--file', '-f', required=True, help='FASTQ file')
    p_qual.set_defaults(func=cmd_quality)

    # entropy
    p_ent = subparsers.add_parser('entropy', help='Shannon entropy')
    p_ent.add_argument('--sequence', '-s', help='DNA sequence')
    p_ent.add_argument('--file', '-f', help='FASTQ file')
    p_ent.add_argument('-k', type=int, default=3, help='K-mer size (default: 3)')
    p_ent.set_defaults(func=cmd_entropy)

    # complexity
    p_comp = subparsers.add_parser('complexity', help='Sequence complexity')
    p_comp.add_argument('--sequence', '-s', help='DNA sequence')
    p_comp.add_argument('--file', '-f', help='FASTQ file')
    p_comp.set_defaults(func=cmd_complexity)

    # analyze
    p_ana = subparsers.add_parser('analyze', help='Full FASTQ analysis')
    p_ana.add_argument('--file', '-f', required=True, help='FASTQ file')
    p_ana.add_argument('--k-values', default='1,2,3', help='Comma-separated k values')
    p_ana.set_defaults(func=cmd_analyze)

    # batch
    p_batch = subparsers.add_parser('batch', help='Batch process from CSV')
    p_batch.add_argument('-i', '--input', required=True, help='Input CSV')
    p_batch.add_argument('-o', '--output', default='results.csv', help='Output CSV')
    p_batch.add_argument('-k', type=int, default=3, help='K-mer size')
    p_batch.set_defaults(func=cmd_batch)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()
