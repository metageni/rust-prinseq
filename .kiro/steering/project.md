# prinseq-rust Project Rules

## Project
Rust reimplementation of PRINSEQ++ — a high-performance FASTQ/FASTA sequence quality control tool.

## Stack
- Rust 2021 edition
- CLI: `clap` (derive)
- Compression: `flate2`
- Deduplication: `bloomfilter`
- Parallelism: `rayon`
- Tests: `tempfile` (dev)

## Conventions
- No `.unwrap()` in production code — use `?` or `.expect("reason")`
- No `unsafe` blocks without a safety comment
- Prefer `&str` over `String` in function parameters
- Keep `main.rs` as thin CLI glue; logic lives in `reads.rs` and other modules
- Run `cargo check` and `cargo test` before committing
- `target/` and `Cargo.lock` are gitignored
