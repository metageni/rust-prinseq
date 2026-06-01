# Changelog

## [1.1.1] - 2026-06-01

### Added
- Bioconda package available: `conda install -c bioconda prinseq-rust` ([PR](https://github.com/bioconda/bioconda-recipes/pulls))
- README updated with Bioconda installation instructions

### Fixed
- `--VERBOSE 2` output now prints exactly 17 lines matching C++ PRINSEQ++ — `trim_to_len` was incorrectly included in the verbose stat list (it is a Rust-only addition and has no C++ counterpart)

## [1.1.0] - 2026-06-01

### Added
- `--trim_to_len N`: trim each read to at most N bases from the 5' end (addresses C++ issue #21)
- `--trim_tail_left` / `--trim_tail_right` now also trim poly-G tails in addition to poly-A/T, supporting 2-colour Illumina sequencing (addresses C++ issue #15)
- Documented order of operations in `--help`: trimmers run first (left→right→tail→qual→to_len), then filters

### Fixed
- **[Issue #9](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/issues/9)** (incorrect mean qual): The C++ `reads.cpp` loop started at `i > 0`, skipping the first base. The Rust port uses `.bytes().map(...).sum()` which correctly includes all bases — this bug never existed in the Rust port.
- **[Issue #10](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/issues/10)** (LANG/locale crash): C++ crashed with `locale::facet::_S_create_c_locale name not valid` when `LANG=C.UTF-8`. Rust does not use `std::locale` and is unaffected — this bug never existed in the Rust port.
- **[Issue #11](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/issues/11)** (dust/entropy score scale): Clarified in `--help` that `--lc_dust` and `--lc_entropy` thresholds are in the range **0–1**, not 0–100 as in old prinseq-lite. A DUST score of 0.5 is roughly equivalent to ~20 in old prinseq-lite.
- **[Issue #12](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/issues/12)** (`-noiupac` removes good reads): The C++ `noiupac` check was case-sensitive and rejected lowercase `acgt`. The Rust port explicitly allows both upper and lowercase `ACGTUNacgtun` — this bug never existed in the Rust port.
- **[Issue #13](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/issues/13)** (seq/qual length mismatch in output): The C++ trimming code could leave seq and qual at different lengths. The Rust port always trims both seq and qual by the same amount — this bug never existed in the Rust port.
- **[Issue #14](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus/issues/14)** (`trim_qual_left` direction): Clarified in `--help` that `--trim_qual_left` trims from the **5' end** (left/start of read). The name "left" refers to the 5' end.

## [1.0.0] - 2026-05-31

### Added
- Full Rust port of PRINSEQ++ released as `prinseq-rust` with feature parity to C++ v1.2
- All filters: `min_len`, `max_len`, `min_gc`, `max_gc`, `min_qual_score`, `min_qual_mean`, `ns_max_n`, `noiupac`, `derep`, `lc_entropy`, `lc_dust`
- All trimmers: `trim_left`, `trim_right`, `trim_tail_left`, `trim_tail_right`, `trim_qual_left`, `trim_qual_right`
- Single-end and paired-end (good/single/bad) read routing
- FASTQ and FASTA input; plain and gzip (`.gz`) for both input and output
- Phred+33 and Phred+64 quality score support (`--phred64`)
- Multi-threading via `--threads` (std::thread + mutex-guarded shared reader)
- Bloom filter deduplication (`--derep`)
- VERBOSE modes 0, 1, 2 matching C++ output format exactly
- Optional-value flags: `--lc_entropy`, `--lc_dust`, `--trim_qual_right`, `--trim_qual_left` default to `0.5`, `0.5`, `20`, `20` when given without a value

### Fixed
- `trim_qual_right` / `trim_qual_left`: C++ had a seq/qual desync bug — `copy_seq.erase(b_win, copy_qual.size())` used the post-erase qual length as the erase count, leaving seq longer than qual. Rust port trims both by the same `step` length. This affects `trim_qual_right` and `trim_qual_rule=gt` cases.

### Compatibility
Output verified to match C++ PRINSEQ++ v1.2 on 59 test cases covering all filters, all trimmers, combinations, paired-end, FASTA input/output, and output formats. 57/59 match exactly; the 2 differences are the intentional `trim_qual_right` bug fix above.

### Performance
- ~1.7–1.8× faster than C++ at all scales (100k–10M reads, single and multi-threaded)
- Benchmarked with a realistic filter set on 150 bp reads; see README for full table
