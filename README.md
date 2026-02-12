# SeqDedupe

**Lightweight sequence deduplication and preprocessing for wet-lab biologists.**

SeqDedupe is an R Shiny application that removes duplicate sequences from FASTA and FASTQ files — both nucleotide and protein — with full transparency. It shows you exactly what was kept, what was removed, and why. Designed for biologists who need clean, auditable sequence data before downstream analysis.

## Why SeqDedupe?

Most deduplication happens via command-line tools (`seqkit rmdup`, `awk`, CD-HIT). These work, but they're opaque — you get an output file with no explanation of what changed.

SeqDedupe gives you:

- A **stepwise pipeline summary** showing exactly how many sequences were removed at each stage
- **Duplicate cluster reports** identifying which sequences matched and which was retained
- A **full audit trail** tagging every input sequence as retained, filtered, or removed
- A **reproducibility report** with parameters, versions, and timestamps for your lab notebook

## Features

- **FASTA and FASTQ support** — reads `.fasta`, `.fa`, `.fna`, `.faa`, `.fastq`, `.fq`; preserves quality scores
- **Nucleotide and protein sequences** — auto-detects sequence type; adapts filters (GC% for DNA, X-count for protein) and labels (bp vs aa) accordingly
- **Multi-file upload** — merge and deduplicate across multiple files
- **Pre-deduplication filtering** — min/max length, GC% range (nucleotide), max N/X count, ambiguous residue threshold
- **Hash-based deduplication** — MD5 hashing for efficient O(n) comparison; deduplicate by sequence, header, or both
- **Before/after visualization** — overlaid length distribution histograms
- **Four downloadable outputs** — deduplicated sequences, audit CSV, cluster CSV, run report
- **No Bioconductor dependency** — pure base R + three CRAN packages

## Installation

```r
# Install dependencies
install.packages(c("shiny", "DT", "digest", "ggplot2"))

# Clone and run
# git clone https://github.com/mbaffour/seqdedupe.git
# cd seqdedupe
shiny::runApp("app.R")
```

### Requirements

- R >= 4.0
- Packages: `shiny`, `DT`, `digest`, `ggplot2`

## Quick Start

1. Launch the app: `shiny::runApp("app.R")`
2. Upload one or more FASTA/FASTQ files (an example file is included in `inst/extdata/`)
3. Set filter thresholds (or leave defaults for deduplication only)
4. Click **Run Pipeline**
5. Review the stepwise summary, duplicate clusters, and audit trail
6. Download your cleaned sequences and reports

### Example

```r
# Try with the included test files
shiny::runApp("app.R")
# Upload inst/extdata/example.fasta (nucleotide)
#   12 sequences: 3 duplicates, 1 short junk, 1 ambiguous
# Or upload inst/extdata/example_protein.fasta (protein)
#   8 sequences: 2 duplicates, 1 ambiguous (X-rich), 1 short peptide
```

## Pipeline

```
Input sequences (nucleotide or protein)
    │
    ├── Length filter (min/max bp or aa)
    ├── GC content filter (nucleotide only)
    ├── N/X count filter (N for DNA, X for protein)
    ├── Ambiguous residue filter (max %)
    │
    └── Hash-based deduplication
            │
            ├── By sequence identity (default)
            ├── By header identity
            └── By both
            
Output: deduplicated sequences + reports
```

## Outputs

| File | Description |
|------|-------------|
| `*_deduped.fasta` | Clean, deduplicated sequences |
| `*_audit.csv` | Every input sequence with status: retained / removed_duplicate / filtered |
| `*_clusters.csv` | Duplicate groups with retained header and removed members |
| `*_seqdedupe_report.txt` | Full reproducibility report: parameters, versions, stepwise counts |

## Use Cases

- **Phage genome curation** — clean redundant sequences from isolation campaigns
- **Amplicon preprocessing** — remove short junk and exact duplicates before analysis
- **Protein database building** — deduplicate proteome FASTA files before annotation or HMM searches
- **Database building** — deduplicate before submitting to repositories
- **Assembly QC** — filter and deduplicate contigs before annotation
- **Teaching** — visual, interactive tool for demonstrating sequence preprocessing concepts

## Limitations

- Exact-match deduplication only (no identity thresholds / clustering)
- In-memory processing — for files >500 MB, use `seqkit rmdup` or CD-HIT
- Single-threaded R — large datasets (>100k sequences) may take a minute

## Roadmap

- [ ] Identity-threshold deduplication (e.g., 99% similarity)
- [ ] Reverse complement awareness
- [ ] R package with CLI interface
- [ ] Shiny Server / Docker deployment
- [ ] Batch processing mode

## Citation

If you use SeqDedupe in your work, please cite:

```
SeqDedupe v0.1.0
https://github.com/mbaffour/seqdedupe
```

## License

[MIT](LICENSE)# SeqDedupe
Shiny app in R for deduplication of FASTA sequences
