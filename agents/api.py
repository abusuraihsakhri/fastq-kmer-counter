"""
FastAPI REST API Server for Fastq Kmer Counter.
"""
from typing import Dict, Any, List
from pathlib import Path
import re
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from .base import AuditLogger, PHIGuard, SecurityException
from .models import SystemTaskPayload, ConsensusDossier
from .supervisor import SystemSupervisor
from kmer_counter import count_kmers, kmer_entropy

supervisor = SystemSupervisor(model_provider="mock")

app = FastAPI(
    title="Fastq Kmer Counter API",
    description="Enterprise Distributed Component Platform (Clinical & Biomedical AI)",
    version="3.0.0-ENTERPRISE",
)


class ChatRequest(BaseModel):
    query: str


class KmerRequest(BaseModel):
    sequence: str = Field(min_length=1, max_length=100000)
    k: int = Field(ge=1, le=100)


@app.post("/api/kmers")
def api_kmers(payload: KmerRequest):
    """Return exact, overlapping forward-strand counts for unambiguous DNA windows."""
    if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", payload.sequence, flags=re.IGNORECASE):
        raise HTTPException(status_code=422, detail="Expected an IUPAC DNA sequence.")
    counts = count_kmers([payload.sequence], payload.k)
    return {
        "k": payload.k,
        "unique_kmers": len(counts),
        "total_kmers": sum(counts.values()),
        "entropy": round(kmer_entropy(counts), 4),
        "counts": counts,
    }


@app.get("/health")
def health():
    return {"status": "HEALTHY", "service": "fastq-kmer-counter", "domain": "Clinical & Biomedical AI", "standard": "CAP / CLSI / ISO Standards", "version": "3.0.0-ENTERPRISE"}


@app.get("/metrics")
def metrics():
    return {
        "dossiers_processed_total": len(supervisor.dossier_registry),
        "audit_blocks_total": len(AuditLogger.get_trail()),
        "system_status": "NOMINAL_OPTIMAL"
    }


@app.post("/api/audit")
def api_audit(payload: SystemTaskPayload):
    try:
        dossier = supervisor.process_task(payload)
    except SecurityException:
        raise HTTPException(status_code=400, detail="Sensitive identifiers are not accepted.")
    return dossier.to_dict()


@app.post("/api/chat")
def api_chat(req: ChatRequest):
    try:
        ans = supervisor.query_supervisory_chat(req.query)
        return {"response": ans}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/audit/logs")
def api_audit_logs():
    return {"audit_trail": AuditLogger.get_trail(), "verified": AuditLogger.verify_integrity()}


# Register after API routes so static frontend paths do not intercept the API.
app.mount("/", StaticFiles(directory=str(Path(__file__).resolve().parent.parent / "web"), html=True), name="web")
