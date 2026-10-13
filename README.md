# FASTQ K-mer Counter

A FASTQ sequence-analysis utility with a dependency-free Python core and a lightweight browser application. It counts overlapping DNA k-mers, summarizes sequencing quality, and estimates distribution entropy.

## Browser application

Open `web/index.html` directly in a modern browser, or serve the repository's `web/` directory with any static web server. Paste four-line FASTQ records or select a local FASTQ file, choose k (1–20), then select **Analyze FASTQ**. The interface reports read/base counts, valid and distinct k-mers, Phred+33 summary quality, Shannon entropy, normalized entropy, and a k-mer frequency spectrum. Export all k-mer counts as CSV.

**Privacy:** Browser analysis runs entirely on the device. The application has no telemetry, uploads, external scripts or server-side processing. A local file is read only after selection. The browser implementation is JavaScript, not Python/Pyodide. It follows the core Python counting rules for forward-strand overlapping k-mers, excluding any windows containing ambiguous bases. Maximum browser input is 8 MiB; for larger files, use the Python CLI (which currently reads input into memory).

GitHub Pages publication is automated by `.github/workflows/pages.yml` once GitHub Pages is configured to deploy from GitHub Actions. A live URL is not listed until deployment can be verified.

## Python CLI

The `kmer_counter.py` library uses only the Python standard library (Python 3.10+). Its FASTQ parser accepts four-line records with IUPAC nucleotide codes and checks headers, separators, quality encoding, record completeness and matching sequence/quality lengths. K-mer counting uses only A, C, G and T windows.

```bash
python cli.py count --sequence ATCGATCG -k 3
python cli.py entropy --sequence ATCGATCG -k 2
python cli.py complexity --sequence ATCGATCG
python cli.py analyze --file reads.fastq --k-values 1,2,3
python cli.py quality --file reads.fastq
python cli.py batch -i sequences.csv -o results.csv -k 3
```

Batch CSV input needs a `sequence` or `seq` column. The CLI limits file inputs to 100 MiB and k to at most 100. The analysis is in-memory, not a streaming, parallel or paired-end pipeline.

## Optional FastAPI server

```bash
python -m pip install "fastapi>=0.111,<1" "uvicorn>=0.30,<1" "pydantic>=2,<3"
uvicorn agents.api:app --host 127.0.0.1 --port 8000
```

Navigate to `http://127.0.0.1:8000/` for the same browser interface. `POST /api/kmers` accepts JSON `{"sequence":"ACGTRGT","k":2}` and returns counts, total/unique k-mers and Shannon entropy. `GET /health` reports availability. The optional pre-existing supervisory demo remains exposed through `/api/audit`, `/api/chat` and `/api/audit/logs`. Its fixed thresholds and deterministic mock response are **not** clinical or scientific validation. The audit ledger is a process-local HMAC chain: it detects modification of retained entries but is not persisted across restarts. The demo endpoints have no authentication and should not be exposed to an untrusted network without additional controls.

To run the API in Docker:

```bash
docker build -t fastq-kmer-counter .
docker run --rm -p 127.0.0.1:8000:8000 fastq-kmer-counter
```

Or use `docker compose up --build`. The optional `AUDIT_SECRET_KEY` environment variable configures the in-memory audit signing key; otherwise a random key is generated for the process.

## Tests

```bash
python -m pip install pytest fastapi uvicorn "pydantic>=2,<3" httpx
python -m pytest -q
node --test web/app.test.cjs
```

GitHub Actions runs Python 3.10–3.12 tests, JavaScript tests and a container build check.

## Structure, compatibility and license

- `kmer_counter.py` — FASTQ parsing, quality statistics, k-mer analysis and CSV processing
- `cli.py` — command-line access to the core and optional supervisory demo
- `web/index.html` and `web/app.js` — standalone offline browser interface (current Chrome, Firefox, Safari, Edge)
- `agents/` — optional FastAPI interface and deterministic supervisory demonstration
- `tests/` and `web/app.test.cjs` — Python and JavaScript regression coverage

MIT license: see [LICENSE](LICENSE).
