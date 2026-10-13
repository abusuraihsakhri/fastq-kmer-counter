"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const {parseFastq, analyze, toCSV} = require("./app.js");

test("parses multiple reads and normalizes nucleotides", () => {
  const reads = parseFastq("@r1\r\nacgt\r\n+\r\nIIII\r\n@r2\nGGTT\n+\n!!!!\n");
  assert.equal(reads.length, 2);
  assert.equal(reads[0].sequence, "ACGT");
});

test("rejects truncated and malformed FASTQ", () => {
  assert.throws(() => parseFastq("@r1\nACGT\n+\nIII\n"), /quality/);
  assert.throws(() => parseFastq("@r1\nACGT\n+\n"), /quality|Incomplete/);
  assert.throws(() => parseFastq("@r1\nACGT\n-\nIIII\n"), /separator/);
  assert.throws(() => parseFastq("@r1\nACGT\n+\nIIII\n@r2"), /Incomplete/);
  assert.throws(() => parseFastq("@r1\nAC?T\n+\nIIII\n"), /invalid DNA/);
});

test("counts overlapping k-mers and excludes all ambiguous bases", () => {
  const reads = parseFastq("@r1\nAAAA\n+\nIIII\n@r2\nACNGT\n+\nIIIII\n");
  const result = analyze(reads, 2);
  assert.equal(result.totalWindows, 7);
  assert.equal(result.totalKmers, 5);
  assert.equal(result.skippedWindows, 2);
  assert.deepEqual(result.ordered[0], {word: "AA", count: 3});
  assert.deepEqual(result.ordered.map(({word}) => word), ["AA", "AC", "GT"]);
});

test("computes entropy, quality and CSV deterministically", () => {
  const reads = parseFastq("@r1\nACGT\n+\nIIII\n");
  const result = analyze(reads, 1);
  assert.equal(result.entropy, 2);
  assert.equal(result.normalizedEntropy, 1);
  assert.equal(result.meanQuality, 40);
  assert.equal(result.percentQ30, 100);
  assert.equal(toCSV(result), "kmer,count\nA,1\nC,1\nG,1\nT,1\n");
});

test("rejects invalid k and handles k larger than the read", () => {
  const reads = parseFastq("@r1\nACG\n+\nIII\n");
  assert.throws(() => analyze(reads, 0), /k-mer/);
  assert.throws(() => analyze(reads, 21), /k-mer/);
  assert.equal(analyze(reads, 4).totalKmers, 0);
});
