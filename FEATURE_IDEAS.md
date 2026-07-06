# Feature Ideas

A backlog of proposed enhancements for SeqDedupe. These are design notes only —
none are implemented yet. Each entry lists the motivation, a sketch of the
approach, and the main risks so a contributor can pick one up.

## 1. Identity-threshold (near-duplicate) deduplication

**Motivation.** SeqDedupe currently removes only exact matches. Real datasets
often contain near-identical sequences (sequencing errors, minor variants) that
exact hashing misses. This is the top item on the README roadmap.

**Sketch.** Add an optional "similarity" mode with a user-set identity threshold
(e.g. 97–100%). For each candidate cluster, compare members with a lightweight
distance (k-mer Jaccard for a fast prefilter, then edit distance via
`utils::adist` on the survivors) and collapse those above the threshold. Keep the
exact-hash path as the default so existing behavior is unchanged.

**Risks.** All-vs-all comparison is O(n^2); needs a k-mer prefilter or length
bucketing to stay usable. Belongs behind a clearly labeled non-default toggle so
users understand results are no longer exact.

## 2. Reverse-complement awareness for nucleotide sequences

**Motivation.** A sequence and its reverse complement represent the same
molecule but hash differently, so true duplicates on opposite strands are kept.
Also on the README roadmap.

**Sketch.** When enabled (nucleotide only), compute a canonical key as
`min(hash(seq), hash(revcomp(seq)))` and deduplicate on that. Reverse complement
is a simple base-map + reverse in base R. Report which orientation was retained
in the cluster CSV.

**Risks.** Only meaningful for DNA/RNA — must be disabled for protein and for
`method = "header"`. Ambiguity codes (N, R, Y, ...) need a defined complement map
or should be passed through unchanged.

## 3. Command-line / batch (non-Shiny) entry point

**Motivation.** The core engine (`parse_fasta`, `apply_filters`,
`deduplicate_seqs`, report generation) is already independent of the UI. Exposing
a scriptable entry point would let users run the same pipeline in headless
batch jobs and pipelines — matching the "R package with CLI interface" and
"Batch processing mode" roadmap items.

**Sketch.** Factor the engine functions into a sourced `R/engine.R` and add a
thin `seqdedupe.R` CLI (via `optparse` or base `commandArgs`) that accepts input
files, filter parameters, and an output directory, then writes the same four
outputs the app produces. The Shiny app would source the shared engine so logic
stays in one place.

**Risks.** Mostly refactoring; the main care is keeping the app and CLI in sync
by sharing one engine file rather than duplicating logic.

## 4. Bundled example datasets surfaced in the UI

**Motivation.** `inst/extdata/example.fasta` and `example_protein.fasta` now ship
with the repo (see the README Quick Start). A "Load example" button would let
first-time users try the full pipeline in one click without hunting for a file.

**Sketch.** Add a small `selectInput`/`actionButton` in the sidebar that loads a
bundled example via `system.file("extdata", ..., package = ...)` (or a relative
path when run from source) and feeds it into the existing pipeline as if it were
uploaded.

**Risks.** Low. Needs a path that resolves both when run from source
(`shiny::runApp`) and when installed as a package.

## 5. Session-persistent run history and report re-download

**Motivation.** Users often run several parameter sets in one sitting and want to
compare them. Today each run replaces the last, so earlier reports are lost.

**Sketch.** Keep a `reactiveVal` list of past runs (parameters + summary counts +
timestamp) within the session and render it as a `DT` table, with a button to
re-download the report/audit for any prior run.

**Risks.** In-memory only (history is lost when the session ends, which is
acceptable). Watch memory use for very large inputs — store summaries and report
text, not full sequence data frames, for old runs.
