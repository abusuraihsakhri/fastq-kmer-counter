"use strict";

/* Lightweight browser implementation of the core Python FASTQ/k-mer workflow. */
const FASTQCounter = (() => {
  const MAX_BYTES = 8 * 1024 * 1024;

  function parseFastq(text) {
    if (typeof text !== "string" || !text.trim()) {
      throw new Error("Paste FASTQ records or choose a FASTQ file.");
    }
    const lines = text.split(/\r?\n/);
    while (lines.length && lines[lines.length - 1] === "") lines.pop();
    if (lines.length % 4 !== 0) {
      throw new Error("Incomplete FASTQ input: each read requires four lines.");
    }
    const records = [];
    for (let i = 0; i < lines.length; i += 4) {
      const [header, sequence, separator, quality] = lines.slice(i, i + 4);
      const number = i / 4 + 1;
      if (!header.startsWith("@") || header.length === 1) {
        throw new Error("Read " + number + ": missing @header.");
      }
      if (!separator.startsWith("+") ||
          (separator.length > 1 && separator.slice(1) !== header.slice(1))) {
        throw new Error("Read " + number + ": invalid + separator.");
      }
      if (!/^[ACGTRYSWKMBDHVN]+$/i.test(sequence)) {
        throw new Error("Read " + number + ": invalid DNA sequence.");
      }
      if (sequence.length !== quality.length || !/^[!-~]+$/.test(quality)) {
        throw new Error("Read " + number + ": quality length or encoding is invalid.");
      }
      records.push({header: header.slice(1), sequence: sequence.toUpperCase(), quality});
    }
    return records;
  }

  function analyze(records, k) {
    if (!Number.isInteger(k) || k < 1 || k > 20) {
      throw new Error("Choose a k-mer length from 1 to 20.");
    }
    const counts = new Map();
    let totalBases = 0, qSum = 0, q20 = 0, q30 = 0, totalWindows = 0;
    for (const record of records) {
      const seq = record.sequence.toUpperCase();
      totalBases += seq.length;
      for (let j = 0; j < record.quality.length; j++) {
        const q = record.quality.charCodeAt(j) - 33;
        qSum += q;
        if (q >= 20) q20++;
        if (q >= 30) q30++;
      }
      for (let j = 0; j <= seq.length - k; j++) {
        const word = seq.slice(j, j + k);
        totalWindows++;
        if (!/^[ACGT]+$/.test(word)) continue;
        counts.set(word, (counts.get(word) || 0) + 1);
      }
    }
    const ordered = Array.from(counts, ([word, count]) => ({word, count}))
      .sort((a, b) => b.count - a.count || a.word.localeCompare(b.word));
    const totalKmers = ordered.reduce((s, x) => s + x.count, 0);
    let entropy = 0;
    const spectrum = new Map();
    for (const {count} of ordered) {
      const p = count / totalKmers;
      entropy -= p * Math.log2(p);
      spectrum.set(count, (spectrum.get(count) || 0) + 1);
    }
    return {
      k, reads: records.length, totalBases, totalWindows, totalKmers,
      skippedWindows: totalWindows - totalKmers, uniqueKmers: counts.size,
      meanQuality: totalBases ? qSum / totalBases : 0,
      percentQ20: totalBases ? 100 * q20 / totalBases : 0,
      percentQ30: totalBases ? 100 * q30 / totalBases : 0,
      entropy, normalizedEntropy: counts.size > 1 ? entropy / Math.log2(counts.size) : 0,
      ordered, spectrum: Array.from(spectrum).sort((a, b) => a[0] - b[0])
    };
  }

  function toCSV(result) {
    return "kmer,count\n" + result.ordered.map(({word, count}) => word + "," + count).join("\n") + "\n";
  }

  function render(result) {
    const stats = document.getElementById("stats");
    const values = [
      ["Reads", result.reads], ["Bases", result.totalBases.toLocaleString()],
      ["Valid k-mers", result.totalKmers.toLocaleString()],
      ["Unique k-mers", result.uniqueKmers.toLocaleString()],
      ["Mean Phred Q", result.meanQuality.toFixed(2)],
      ["Shannon entropy (bits)", result.entropy.toFixed(4)],
      ["Normalized entropy", result.normalizedEntropy.toFixed(4)],
      ["Bases ≥ Q30", result.percentQ30.toFixed(1) + "%"]
    ];
    stats.replaceChildren();
    for (const [name, value] of values) {
      const item = document.createElement("div");
      item.className = "stat";
      const dt = document.createElement("span"), dd = document.createElement("strong");
      dt.textContent = name; dd.textContent = value;
      item.append(dt, dd); stats.append(item);
    }
    const body = document.getElementById("counts");
    body.replaceChildren();
    for (const {word, count} of result.ordered.slice(0, 25)) {
      const row = document.createElement("tr");
      const wordCell = document.createElement("td"), countCell = document.createElement("td");
      wordCell.textContent = word; countCell.textContent = count.toLocaleString();
      row.append(wordCell, countCell); body.append(row);
    }
    const spectrum = document.getElementById("spectrum");
    spectrum.replaceChildren();
    for (const [frequency, n] of result.spectrum.slice(0, 30)) {
      const item = document.createElement("div");
      item.className = "spectrum-item";
      item.textContent = frequency + "×: " + n + " unique";
      spectrum.append(item);
    }
    document.getElementById("results").hidden = false;
    document.getElementById("download").disabled = false;
    document.getElementById("notes").textContent =
      result.skippedWindows.toLocaleString() + " windows with ambiguous bases excluded. " +
      "Showing the 25 most frequent k-mers; CSV includes all " + result.uniqueKmers.toLocaleString() + ".";
  }

  function init() {
    const form = document.getElementById("counterForm");
    const upload = document.getElementById("upload");
    const input = document.getElementById("fastq");
    const status = document.getElementById("status");
    let lastResult = null;
    form.addEventListener("submit", async (event) => {
      event.preventDefault();
      const run = document.getElementById("run");
      run.disabled = true;
      status.className = "";
      status.textContent = "Analyzing reads…";
      document.getElementById("results").hidden = true;
      document.getElementById("download").disabled = true;
      try {
        const file = upload.files && upload.files[0];
        if (file && file.size > MAX_BYTES) {
          throw new Error("File exceeds the 8 MiB browser limit; use the Python CLI for larger files.");
        }
        const text = file ? await file.text() : input.value;
        if (new TextEncoder().encode(text).length > MAX_BYTES) {
          throw new Error("Input exceeds the 8 MiB browser limit.");
        }
        const k = Number(document.getElementById("k").value);
        const result = analyze(parseFastq(text), k);
        lastResult = result;
        render(result);
        status.textContent = "Analysis complete. Processing stayed in this browser.";
        status.className = "ok";
      } catch (error) {
        lastResult = null;
        status.textContent = error.message || "Unable to analyze this FASTQ input.";
        status.className = "error";
      } finally {
        run.disabled = false;
      }
    });
    document.getElementById("download").addEventListener("click", () => {
      if (!lastResult) return;
      const blob = new Blob([toCSV(lastResult)], {type: "text/csv;charset=utf-8"});
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url; link.download = "kmers-k" + lastResult.k + ".csv";
      document.body.append(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    });
    upload.addEventListener("change", () => {
      const name = upload.files && upload.files[0] ? upload.files[0].name : "Pasted FASTQ text";
      document.getElementById("source").textContent = name;
    });
    document.getElementById("clearFile").addEventListener("click", () => {
      upload.value = "";
      document.getElementById("source").textContent = "Pasted FASTQ text";
    });
  }

  return {MAX_BYTES, parseFastq, analyze, toCSV, init};
})();

if (typeof module !== "undefined" && module.exports) module.exports = FASTQCounter;
if (typeof document !== "undefined") document.addEventListener("DOMContentLoaded", FASTQCounter.init);
