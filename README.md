# Fastq Kmer Counter

> **Domain:** Bioinformatics & Genomic Sequence Analysis

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

FASTQ K-mer Counter parses FASTQ files, counts k-mers for arbitrary k, computes entropy,
quality score statistics (Phred), and sequence complexity estimation.

The core library (`kmer_counter.py`) has **zero external dependencies** — Python stdlib only.
Optional FastAPI server and agent framework require `fastapi`, `uvicorn`, and `pydantic`.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`parse_fastq()`**: Parse FASTQ format text into list of {header, sequence, quality} dicts.
- **`parse_fastq_file()`**: Parse a FASTQ file from disk (with path traversal protection).
- **`count_kmers()`**: Count all k-mers across a list of sequences.
- **`count_kmers_single()`**: Count k-mers in a single sequence.
- **`kmer_frequency_spectrum()`**: Compute frequency spectrum: histogram of k-mer counts.
- **`kmer_entropy()`**: Shannon entropy of k-mer distribution.
- **`normalized_entropy()`**: Normalized entropy (0-1 scale).
- **`per_read_quality()`**: Quality statistics for a single read.
- **`per_position_quality()`**: Mean quality at each position across all reads.
- **`overall_quality_stats()`**: Overall quality statistics across all reads.
- **`sequence_complexity()`**: Estimate sequence complexity using entropy, dinucleotide diversity, and homopolymer detection.
- **`analyze_fastq()`**: Complete analysis combining all metrics.
- **`process_batch()`**: Batch process sequences from CSV.

---

## 💻 CLI Quickstart & Usage

### Installation
```bash
# Core functionality requires no dependencies
# For optional API server:
pip install fastapi uvicorn pydantic pytest
```

### 1. Count K-mers
```bash
python cli.py count --sequence ATCGATCG -k 3
python cli.py count --file sample.fastq -k 3 --top 10
```

### 2. Quality Statistics
```bash
python cli.py quality --file sample.fastq
```

### 3. Entropy Analysis
```bash
python cli.py entropy --sequence ATCGATCG -k 2
python cli.py entropy --file sample.fastq -k 3
```

### 4. Sequence Complexity
```bash
python cli.py complexity --sequence ATCGATCGATCG
python cli.py complexity --file sample.fastq
```

### 5. Full Analysis
```bash
python cli.py analyze --file sample.fastq --k-values 1,2,3
```

### 6. Batch Processing
```bash
python cli.py batch -i input.csv -o output.csv -k 3
```

### 7. Audit & Verification (Agent Framework)
```bash
python cli.py audit --task-id TASK-001 --primary-metric 14.5
python cli.py chat "Explain k-mer counting"
python cli.py verify-audit
```

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation.
* **Path Traversal Protection:** All file operations validate paths against the project directory.
* **Input Validation:** K-mer parameters validated to prevent resource exhaustion.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Run specific test files:
```bash
pytest test_kmer_counter.py -v        # Core functionality (36 tests)
pytest tests/test_enrichment.py -v    # Enrichment modules (2 tests)
pytest tests/test_fastq_kmer_counter.py -v  # Agent framework (3 tests)
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t fastq-kmer-counter .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret-key fastq-kmer-counter
```

Or with docker-compose:
```bash
AUDIT_SECRET_KEY=your-secret-key docker-compose up
```

---

## 📁 Project Structure

```
fastq-kmer-counter/
├── kmer_counter.py      # Core FASTQ analysis library
├── cli.py               # Command-line interface
├── enrichment.py        # Enrichment feature engines
├── simulator.py         # High-throughput simulation
├── agents/              # Agent framework
│   ├── base.py          # PHI guard, HMAC audit trail
│   ├── models.py        # Pydantic data models
│   ├── supervisor.py    # Supervisor orchestrator
│   ├── workers.py       # Specialized domain workers
│   ├── api.py           # FastAPI REST endpoints
│   ├── metrics.py       # Prometheus metrics
│   ├── llm_factory.py   # LLM provider factory
│   ├── learning.py      # Bayesian calibration engine
│   └── streamer.py      # WebSocket telemetry
├── tests/               # Test suite
├── web/                 # Operations console (HTML)
├── Dockerfile           # Container definition
├── docker-compose.yml   # Compose configuration
└── openapi_spec.json    # OpenAPI specification
```
