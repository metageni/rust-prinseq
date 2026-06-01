# rust-prinseq

A Rust port of [PRINSEQ++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) — a fast FASTQ/FASTA quality control, filtering, and trimming tool for genomic and metagenomic sequence data.

## Performance

Benchmarked against the original C++ on 150 bp reads with a realistic filter set:
`--min_len 50 --max_len 200 --min_qual_mean 20 --ns_max_n 5 --lc_entropy=0.5 --trim_qual_right=20 --trim_tail_right=5`

| Input size | C++ 1 thread | Rust 1 thread | Speedup | C++ 4 threads | Rust 4 threads | Speedup |
|---|---|---|---|---|---|---|
| 100k reads (30 MB) | 0.79s | 0.47s | 1.7× | 0.30s | 0.19s | 1.6× |
| 500k reads (150 MB) | 3.93s | 2.35s | 1.7× | 1.46s | 0.83s | 1.8× |
| 1M reads (301 MB) | 7.91s | 4.64s | 1.7× | 2.84s | 1.62s | 1.8× |
| 10M reads (3 GB) | 78.79s | 46.78s | 1.7× | 29.05s | 16.58s | 1.8× |

**rust-prinseq is consistently ~1.7–1.8× faster** than C++ PRINSEQ++ at all scales, both single and multi-threaded.
Multi-threading scales similarly in both tools (~1.7–2× gain from 1→4 threads).

## Installation

### Prerequisites

You need the Rust toolchain (version 1.70 or later). If you don't have it:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env   # or restart your shell
```

Verify:

```bash
rustc --version   # should print rustc 1.70.0 or later
cargo --version
```

### Option 1 — Install from source (recommended)

```bash
git clone https://github.com/metageni/rust-prinseq.git
cd rust-prinseq
cargo install --path .
```

This compiles a release binary and places it in `~/.cargo/bin/`, which is already on your `PATH` after the Rust install. You can then call `rust-prinseq` from anywhere.

### Option 2 — Build manually

```bash
git clone https://github.com/metageni/rust-prinseq.git
cd rust-prinseq
cargo build --release
```

The binary is at `target/release/rust-prinseq`. Copy it wherever you like:

```bash
cp target/release/rust-prinseq /usr/local/bin/   # system-wide
# or
cp target/release/rust-prinseq ~/bin/             # user-local
```

### Verify the installation

```bash
rust-prinseq --version
rust-prinseq --help
```

### Uninstall

If installed via `cargo install`:

```bash
cargo uninstall rust-prinseq
```

## Usage

```
rust-prinseq --fastq <file> [OPTIONS]
```

### Input

| Flag | Description |
|---|---|
| `--fastq <file>` | Input FASTQ (plain or `.gz`) |
| `--fastq2 <file>` | Second FASTQ for paired-end |
| `--FASTA` | Input is FASTA (quality set to A for all bases) |
| `--phred64` | Input quality is Phred+64 (old Illumina/Solexa) |

### Output

| Flag | Description |
|---|---|
| `--out_name <str>` | Base name for output files (default: random 6-char string) |
| `--out_format <0\|1>` | 0 = FASTQ (default), 1 = FASTA |
| `--out_gz` | Write gzip-compressed output |
| `--rm_header` | Strip header from the `+` line |
| `--out_good`, `--out_bad`, `--out_single` | Override individual output file paths (R1) |
| `--out_good2`, `--out_bad2`, `--out_single2` | Override individual output file paths (R2) |

For single-end, output files are `<out_name>_good_out.fastq` and `<out_name>_bad_out.fastq`.  
For paired-end, output files are `_good_out_R1`, `_good_out_R2`, `_single_out_R1`, `_single_out_R2`, `_bad_out_R1`, `_bad_out_R2`.

### Filters

| Flag | Description |
|---|---|
| `--min_len <int>` | Remove reads shorter than N |
| `--max_len <int>` | Remove reads longer than N |
| `--min_gc <float>` | Remove reads with GC% below N |
| `--max_gc <float>` | Remove reads with GC% above N |
| `--min_qual_score <int>` | Remove reads with any base quality below N |
| `--min_qual_mean <int>` | Remove reads with mean quality below N |
| `--ns_max_n <int>` | Remove reads with more than N Ns |
| `--noiupac` | Remove reads with non-ACGTUN characters |
| `--derep` | Remove exact duplicate sequences (bloom filter) |
| `--lc_entropy[=float]` | Remove low-complexity reads by entropy (default threshold: 0.5) |
| `--lc_dust[=float]` | Remove low-complexity reads by DUST score (default threshold: 0.5) |

### Trimming

| Flag | Description |
|---|---|
| `--trim_left <int>` | Trim N bases from 5' end |
| `--trim_right <int>` | Trim N bases from 3' end |
| `--trim_tail_left <int>` | Trim poly-A/T tail from 5' end (min length N) |
| `--trim_tail_right <int>` | Trim poly-A/T tail from 3' end (min length N) |
| `--trim_qual_left[=float]` | Quality-trim 5' end (default threshold: 20) |
| `--trim_qual_right[=float]` | Quality-trim 3' end (default threshold: 20) |
| `--trim_qual_type <str>` | Score type: `min`, `mean` (default), `max`, `sum` |
| `--trim_qual_rule <str>` | Rule: `lt` (default), `gt`, `et` |
| `--trim_qual_window <int>` | Window size (default: 5) |
| `--trim_qual_step <int>` | Step size (default: 2) |

### Other

| Flag | Description |
|---|---|
| `--threads <int>` | Number of threads (default: 1) |
| `--VERBOSE <0\|1\|2>` | 0 = silent, 1 = active filters only (default), 2 = all counts |

## Examples

```bash
# Basic quality filter
rust-prinseq --fastq reads.fastq --min_len 50 --min_qual_mean 20 --out_name clean

# Paired-end with trimming
rust-prinseq --fastq R1.fastq.gz --fastq2 R2.fastq.gz \
  --min_len 50 --trim_qual_right=20 --trim_tail_right=5 \
  --threads 4 --out_name clean

# Low-complexity filter + deduplication
rust-prinseq --fastq reads.fastq --lc_entropy --derep --out_name clean

# FASTA input, FASTA output
rust-prinseq --fastq reads.fasta --FASTA --out_format 1 --min_len 100 --out_name clean
```

## Differences from C++ PRINSEQ++

- No autoconf/Makefile — standard `cargo build`
- No Boost dependency
- `--lc_entropy` / `--lc_dust` / `--trim_qual_right` / `--trim_qual_left` accept an optional `=value` (e.g. `--lc_entropy=0.7`) or no value (uses default)
- Fixes a minor bug in C++ `trim_qual_right` where sequence and quality could desync after trimming

## Citation

If you use this tool, please also cite the original PRINSEQ++ paper:

> Cantu VA, Sadural J, Edwards R. *PRINSEQ++, a multi-threaded tool for fast and efficient quality control and preprocessing of sequencing datasets.* PeerJ Preprints, 2019.
