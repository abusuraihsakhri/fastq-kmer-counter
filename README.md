# Fastq Kmer Counter

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

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

FASTQ K-mer Counter — Real Implementation

Parses FASTQ files, counts k-mers for arbitrary k, computes entropy,
quality score statistics (Phred), and sequence complexity estimation.

Zero external dependencies — Python stdlib only.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`parse_fastq()`**: Parse FASTQ format text into list of {header, sequence, quality} dicts.

FASTQ format: 4 lines per record:
    Line 1: @header
    Line 2: sequence
    Line 3: + (optional repeated header)
    Line 4: quality string (Phred+33 encoded)
- **`parse_fastq_file()`**: Parse a FASTQ file from disk.
- **`count_kmers()`**: Count all k-mers across a list of sequences.

Args:
    sequences: List of DNA sequences.
    k: K-mer length.

Returns:
    Dict mapping k-mer string to count.
- **`count_kmers_single()`**: Count k-mers in a single sequence.
- **`kmer_frequency_spectrum()`**: Compute frequency spectrum: histogram of k-mer counts.

Returns dict mapping count_value -> number_of_kmers_with_that_count.

---

## 📐 Mathematical Formulation & Logic

```text
  """Calculate Shannon entropy of k-mer distribution.
  """Calculate quality statistics for a single read."""
  """Calculate mean quality at each position across all reads.
  """Calculate overall quality statistics across all reads."""
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --file <value> --sequence <value> --top <value> --k-values <value>
```

### Parameter Reference
- `--file`: Specifies input measurement or parameter value.
- `--sequence`: Specifies input measurement or parameter value.
- `--top`: Specifies input measurement or parameter value.
- `--k-values`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `query` | Parameter / observation metric | Required |
| `name` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t fastq-kmer-counter .
docker run -p 8000:8000 fastq-kmer-counter
```
