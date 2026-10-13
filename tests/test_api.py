"""HTTP smoke tests for the optional FastAPI server and static frontend."""
from fastapi.testclient import TestClient
from agents.api import app

client = TestClient(app)


def test_health_and_frontend():
    assert client.get("/health").status_code == 200
    page = client.get("/")
    assert page.status_code == 200
    assert "FASTQ K-mer Counter" in page.text
    script = client.get("/app.js")
    assert script.status_code == 200
    assert "parseFastq" in script.text


def test_api_counts_and_ambiguity_handling():
    response = client.post("/api/kmers", json={"sequence": "ACRGT", "k": 2})
    assert response.status_code == 200
    assert response.json()["counts"] == {"AC": 1, "GT": 1}
    assert response.json()["total_kmers"] == 2


def test_api_input_validation():
    assert client.post("/api/kmers", json={"sequence": "ACGT", "k": 0}).status_code == 422
    assert client.post("/api/kmers", json={"sequence": "not dna", "k": 2}).status_code == 422


def test_api_blocks_identifying_information():
    response = client.post("/api/audit", json={
        "task_id": "MRN-123456", "target_identifier": "SAMPLE",
        "primary_metric": 4.0,
    })
    assert response.status_code == 400
    assert "Sensitive identifiers" in response.text
