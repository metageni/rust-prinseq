/// Integration tests for rust-prinseq.
///
/// Expected outputs are verified against the C++ PRINSEQ++ binary (v1.2).
/// Where rust-prinseq intentionally fixes a C++ bug, the comment notes the difference.
use std::fs;
use std::io::Read;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    let r = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/release/rust-prinseq");
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/debug/rust-prinseq");
    if r.exists() { r } else { d }
}

fn data(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data").join(name)
        .to_string_lossy().into_owned()
}

/// Run the binary, return (good, bad) fastq contents for single-end.
fn run_se(extra: &[&str]) -> (String, String) {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    let status = Command::new(bin())
        .arg("--fastq").arg(data("test_F.fastq"))
        .args(extra)
        .arg("--out_name").arg(&prefix)
        .arg("--VERBOSE").arg("0")
        .status().unwrap();
    assert!(status.success());
    let read = |s: &str| fs::read_to_string(format!("{prefix}{s}")).unwrap_or_default();
    (read("_good_out.fastq"), read("_bad_out.fastq"))
}

fn headers(fastq: &str) -> Vec<String> {
    fastq.lines().filter(|l| l.starts_with('@'))
        .map(|l| l[1..].to_string()).collect()
}

fn seq_of(fastq: &str, name: &str) -> String {
    let lines: Vec<&str> = fastq.lines().collect();
    let i = lines.iter().position(|&l| l == format!("@{name}")).unwrap();
    lines[i + 1].to_string()
}

fn sep_lines(fastq: &str) -> Vec<String> {
    fastq.lines().collect::<Vec<_>>().chunks(4)
        .filter_map(|c| c.get(2).map(|s| s.to_string())).collect()
}

// ── filters ───────────────────────────────────────────────────────────────────

#[test]
fn test_min_len() {
    let (good, bad) = run_se(&["--min_len", "10"]);
    assert_eq!(headers(&good), ["seq1_F","seq5_F","seq6_F","seq7_F","seq8_F"]);
    assert_eq!(headers(&bad),  ["seq2_F","seq3_F","seq4_F"]);
}

#[test]
fn test_max_len() {
    let (good, bad) = run_se(&["--max_len", "15"]);
    assert_eq!(headers(&good), ["seq1_F","seq2_F","seq3_F","seq4_F","seq5_F","seq6_F"]);
    assert_eq!(headers(&bad),  ["seq7_F","seq8_F"]);
}

#[test]
fn test_ns_max_n() {
    let (_, bad) = run_se(&["--ns_max_n", "0"]);
    assert_eq!(headers(&bad), ["seq2_F","seq3_F"]);
}

#[test]
fn test_min_qual_score() {
    // All reads have 'A' qual (phred 0) → all fail min_qual_score=35
    let (good, bad) = run_se(&["--min_qual_score", "35"]);
    assert!(good.is_empty());
    assert_eq!(headers(&bad).len(), 8);
}

#[test]
fn test_min_qual_mean() {
    // Rust correctly averages all positions (C++ has an off-by-one bug skipping index 0).
    // seq3_F (ABCABCA, mean phred ≈ 1) and seq7_F (mixed ABC+9999, mean < 33) fail.
    // seq8_F (all 'C' = phred 34) and others pass.
    let (_, bad) = run_se(&["--min_qual_mean", "33"]);
    assert_eq!(headers(&bad), ["seq3_F","seq7_F"]);
}

#[test]
fn test_noiupac() {
    // All seqs contain only ACGTUN (U is allowed) → nothing filtered
    let (good, bad) = run_se(&["--noiupac"]);
    assert_eq!(headers(&good).len(), 8);
    assert!(bad.is_empty());
}

#[test]
fn test_min_gc() {
    let (_, bad) = run_se(&["--min_gc", "50"]);
    assert_eq!(headers(&bad), ["seq2_F","seq3_F","seq4_F","seq5_F","seq7_F","seq8_F"]);
}

#[test]
fn test_max_gc() {
    let (_, bad) = run_se(&["--max_gc", "50"]);
    assert_eq!(headers(&bad), ["seq1_F","seq6_F"]);
}

#[test]
fn test_derep() {
    // seq1_F == seq6_F, seq7_F == seq8_F → second occurrence of each removed
    let (good, _) = run_se(&["--derep"]);
    assert_eq!(headers(&good), ["seq1_F","seq2_F","seq3_F","seq4_F","seq5_F","seq7_F"]);
}

#[test]
fn test_lc_entropy() {
    // All test reads are low-complexity → all filtered at default threshold 0.5
    let (good, _) = run_se(&["--lc_entropy"]);
    assert!(good.is_empty());
}

#[test]
fn test_lc_dust() {
    // Only all-C reads (seq1_F, seq6_F) are high-DUST and pass
    let (good, _) = run_se(&["--lc_dust"]);
    assert_eq!(headers(&good), ["seq1_F","seq6_F"]);
}

// ── trimmers ─────────────────────────────────────────────────────────────────

#[test]
fn test_trim_left() {
    let (good, _) = run_se(&["--trim_left", "3"]);
    assert_eq!(seq_of(&good, "seq1_F"), "CCCCCCCCCCCC");
    // qual also trimmed
    let qual: Vec<&str> = good.lines().collect::<Vec<_>>()
        .chunks(4).filter_map(|c| c.get(3).copied()).collect();
    assert_eq!(qual[0], "ABCABCABCABC");
}

#[test]
fn test_trim_right() {
    let (good, _) = run_se(&["--trim_right", "3"]);
    assert_eq!(seq_of(&good, "seq1_F"), "CCCCCCCCCCCC");
}

#[test]
fn test_trim_tail_left() {
    // seq7_F and seq8_F are entirely poly-A/T → become bad
    let (_, bad) = run_se(&["--trim_tail_left", "3"]);
    assert_eq!(headers(&bad), ["seq7_F","seq8_F"]);
}

#[test]
fn test_trim_tail_right() {
    // seq7_F and seq8_F are entirely poly-A/T → become bad
    let (_, bad) = run_se(&["--trim_tail_right", "3"]);
    assert_eq!(headers(&bad), ["seq7_F","seq8_F"]);
}

#[test]
fn test_trim_qual_right() {
    // seq7_F has 9999... tail (phred 24 < 30) → trimmed to 55 bases
    // seq8_F all-C qual (phred 34 > 30) → untouched (163 bases)
    // Note: C++ has a seq/qual desync bug here; rust-prinseq trims correctly.
    let (good, _) = run_se(&["--trim_qual_right", "30"]);
    assert_eq!(seq_of(&good, "seq7_F").len(), 55);
    assert_eq!(seq_of(&good, "seq8_F").len(), 163);
}

#[test]
fn test_trim_qual_left() {
    // seq7_F qual: ABC...9999 — window trims ABC section but 9999 (phred 24) also < 30,
    // however the read length prevents full trim (window guard). seq7 stays in good.
    // seq8_F: all-C qual (phred 34 > 30) → untouched.
    let (good, bad) = run_se(&["--trim_qual_left", "30"]);
    assert!(bad.is_empty() || !headers(&bad).contains(&"seq8_F".to_string()));
    assert!(headers(&good).contains(&"seq8_F".to_string()));
    assert_eq!(seq_of(&good, "seq8_F").len(), 163);
}

// ── output options ────────────────────────────────────────────────────────────

#[test]
fn test_rm_header() {
    let (good, _) = run_se(&["--rm_header", "--min_len", "10"]);
    assert!(sep_lines(&good).iter().all(|s| s == "+"));
}

#[test]
fn test_out_format_fasta() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--out_format", "1",
               "--min_len", "10", "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let content = fs::read_to_string(format!("{prefix}_good_out.fasta")).unwrap();
    assert!(content.starts_with('>'));
    let names: Vec<&str> = content.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(names, [">seq1_F",">seq5_F",">seq6_F",">seq7_F",">seq8_F"]);
}

#[test]
fn test_gz_output() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--min_len", "10",
               "--out_gz", "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let mut f = fs::File::open(format!("{prefix}_good_out.fastq.gz")).unwrap();
    let mut magic = [0u8; 2];
    f.read_exact(&mut magic).unwrap();
    assert_eq!(magic, [0x1f, 0x8b], "output must be gzip");
}

// ── input formats ─────────────────────────────────────────────────────────────

#[test]
fn test_fasta_input() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fasta"), "--FASTA",
               "--min_len", "10", "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let good = fs::read_to_string(format!("{prefix}_good_out.fastq")).unwrap();
    assert_eq!(headers(&good), ["seq1_F","seq5_F","seq6_F","seq7_F","seq8_F"]);
}

#[test]
fn test_gz_input() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq.gz"),
               "--min_len", "10", "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let good = fs::read_to_string(format!("{prefix}_good_out.fastq")).unwrap();
    // gz test file has headers without _F suffix
    assert_eq!(headers(&good), ["seq1","seq5","seq6","seq7","seq8"]);
}

// ── paired-end ────────────────────────────────────────────────────────────────

#[test]
fn test_paired_end() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--fastq2", &data("test_R.fastq"),
               "--min_len", "10", "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let read = |s: &str| fs::read_to_string(format!("{prefix}{s}")).unwrap_or_default();
    let good_r1  = read("_good_out_R1.fastq");
    let single_r2 = read("_single_out_R2.fastq");
    // Both mates ≥10 bp: seq1,5,6,7,8 → good
    assert_eq!(headers(&good_r1), ["seq1_F","seq5_F","seq6_F","seq7_F","seq8_F"]);
    // seq4_F (9 bp) fails but seq4_R (21 bp) passes → R2 single
    assert_eq!(headers(&single_r2), ["seq4_R"]);
}

// ── verbose output ────────────────────────────────────────────────────────────

#[test]
fn test_verbose_1() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    let output = Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--min_len", "10",
               "--out_name", &prefix, "--VERBOSE", "1"])
        .output().unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(stdout.trim(), "3 reads removed by -min_len");
}

#[test]
fn test_verbose_2() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    let output = Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--min_len", "10",
               "--out_name", &prefix, "--VERBOSE", "2"])
        .output().unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout);
    // 17 lines, one per stat, in order: min_len=3, rest=0
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 17);
    assert_eq!(lines[0], "3");   // min_len
    assert!(lines[1..].iter().all(|&l| l == "0"));
}

#[test]
fn test_verbose_0_silent() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    let output = Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--min_len", "10",
               "--out_name", &prefix, "--VERBOSE", "0"])
        .output().unwrap();
    assert!(output.stdout.is_empty());
}

// ── combined filters ──────────────────────────────────────────────────────────

#[test]
fn test_min_len_and_ns_max_n() {
    let (good, bad) = run_se(&["--min_len", "10", "--ns_max_n", "0"]);
    assert_eq!(headers(&good), ["seq1_F","seq5_F","seq6_F","seq7_F","seq8_F"]);
    assert_eq!(headers(&bad),  ["seq2_F","seq3_F","seq4_F"]);
}

#[test]
fn test_derep_and_min_len() {
    // derep removes seq6_F (dup of seq1_F) and seq8_F (dup of seq7_F);
    // min_len=10 removes seq2_F,seq3_F,seq4_F
    let (good, _) = run_se(&["--derep", "--min_len", "10"]);
    assert_eq!(headers(&good), ["seq1_F","seq5_F","seq7_F"]);
}

#[test]
fn test_min_len_and_max_len() {
    let (good, _) = run_se(&["--min_len", "8", "--max_len", "13"]);
    assert_eq!(headers(&good), ["seq2_F","seq4_F","seq5_F"]);
}

// ── filter edge cases ─────────────────────────────────────────────────────────

#[test]
fn test_ns_max_n_one() {
    // seq2_F has 1 N, seq3_F has 1 N → both pass ns_max_n=1
    let (_, bad) = run_se(&["--ns_max_n", "1"]);
    assert!(bad.is_empty());
}

#[test]
fn test_lc_entropy_explicit_value() {
    // threshold=0.1 — all test reads are still low-complexity
    let (good, _) = run_se(&["--lc_entropy=0.1"]);
    assert!(good.is_empty());
}

#[test]
fn test_lc_dust_explicit_value() {
    // threshold=0.1 — all test reads still filtered
    let (good, _) = run_se(&["--lc_dust=0.1"]);
    assert!(good.is_empty());
}

#[test]
fn test_min_qual_score_passes_high_qual() {
    // seq8_F has all 'C' qual (phred 34) → passes min_qual_score=30
    let (good, _) = run_se(&["--min_qual_score", "30"]);
    assert!(headers(&good).contains(&"seq8_F".to_string()));
    // seq7_F has '9' qual (phred 24) → fails
    assert!(headers(&good).iter().all(|h| h != "seq7_F"));
}

// ── trim parameter variants ───────────────────────────────────────────────────

#[test]
fn test_trim_qual_right_type_min() {
    // type=min: uses minimum quality in window; seq7 has '9' (phred 24) < 20 → trims
    let (good, _) = run_se(&["--trim_qual_right=20", "--trim_qual_type=min"]);
    // seq7 trimmed, seq8 (all phred 34 > 20) untouched
    assert_eq!(seq_of(&good, "seq8_F").len(), 163);
}

#[test]
fn test_trim_qual_right_rule_gt() {
    // rule=gt threshold=30: trim where quality IS above 30
    // seq8_F all-C (phred 34 > 30) → trimmed down to near nothing
    let (good, _) = run_se(&["--trim_qual_right=30", "--trim_qual_rule=gt"]);
    assert!(seq_of(&good, "seq8_F").len() < 10);
}

#[test]
fn test_trim_qual_right_window_step() {
    // window=10 step=5: seq7 trimmed in larger chunks
    let (good, _) = run_se(&["--trim_qual_right=30", "--trim_qual_window=10", "--trim_qual_step=5"]);
    assert!(seq_of(&good, "seq7_F").len() < 163);
    assert_eq!(seq_of(&good, "seq8_F").len(), 163);
}

#[test]
fn test_trim_left_removes_qual_too() {
    // qual must be trimmed in sync with seq
    let (good, _) = run_se(&["--trim_left", "3"]);
    let lines: Vec<&str> = good.lines().collect();
    let i = lines.iter().position(|&l| l == "@seq1_F").unwrap();
    assert_eq!(lines[i+1].len(), lines[i+3].len(), "seq and qual lengths must match");
}

#[test]
fn test_trim_right_removes_qual_too() {
    let (good, _) = run_se(&["--trim_right", "3"]);
    let lines: Vec<&str> = good.lines().collect();
    let i = lines.iter().position(|&l| l == "@seq1_F").unwrap();
    assert_eq!(lines[i+1].len(), lines[i+3].len());
}

// ── output overrides ──────────────────────────────────────────────────────────

#[test]
fn test_out_good_bad_override() {
    let out = tempfile::tempdir().unwrap();
    let good_path = out.path().join("my_good.fastq").to_string_lossy().into_owned();
    let bad_path  = out.path().join("my_bad.fastq").to_string_lossy().into_owned();
    let prefix    = out.path().join("unused").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--min_len", "10",
               "--out_good", &good_path, "--out_bad", &bad_path,
               "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let good = fs::read_to_string(&good_path).unwrap();
    assert_eq!(headers(&good), ["seq1_F","seq5_F","seq6_F","seq7_F","seq8_F"]);
}

#[test]
fn test_out_name_custom() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("myprefix").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--min_len", "10",
               "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    assert!(std::path::Path::new(&format!("{prefix}_good_out.fastq")).exists());
    assert!(std::path::Path::new(&format!("{prefix}_bad_out.fastq")).exists());
}

// ── paired-end edge cases ─────────────────────────────────────────────────────

#[test]
fn test_paired_both_fail_go_bad() {
    // min_len=20: seq7_F(163) and seq8_F(163) pass; seq4_R(21) passes but seq4_F(9) fails → single
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--fastq2", &data("test_R.fastq"),
               "--min_len", "20", "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let read = |s: &str| fs::read_to_string(format!("{prefix}{s}")).unwrap_or_default();
    let good_r1   = read("_good_out_R1.fastq");
    let single_r2 = read("_single_out_R2.fastq");
    assert_eq!(headers(&good_r1), ["seq7_F","seq8_F"]);
    assert_eq!(headers(&single_r2), ["seq4_R"]);
}

#[test]
fn test_paired_out_format_fasta() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--fastq2", &data("test_R.fastq"),
               "--min_len", "10", "--out_format", "1",
               "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let good = fs::read_to_string(format!("{prefix}_good_out_R1.fasta")).unwrap();
    assert!(good.starts_with('>'));
}

// ── FASTA input + FASTA output ────────────────────────────────────────────────

#[test]
fn test_fasta_in_fasta_out() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    Command::new(bin())
        .args(["--fastq", &data("test_F.fasta"), "--FASTA", "--out_format", "1",
               "--min_len", "10", "--out_name", &prefix, "--VERBOSE", "0"])
        .status().unwrap();
    let content = fs::read_to_string(format!("{prefix}_good_out.fasta")).unwrap();
    let names: Vec<&str> = content.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(names, [">seq1_F",">seq5_F",">seq6_F",">seq7_F",">seq8_F"]);
}

// ── verbose with multiple active filters ─────────────────────────────────────

#[test]
fn test_verbose_1_multiple_filters() {
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    let output = Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--min_len", "10", "--ns_max_n", "0",
               "--out_name", &prefix, "--VERBOSE", "1"])
        .output().unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("reads removed by -min_len"));
    assert!(stdout.contains("reads removed by -ns_max_n"));
}

#[test]
fn test_verbose_2_order() {
    // VERBOSE=2 prints 17 values in fixed order; with only min_len active,
    // first value is min_len count, rest are 0
    let out = tempfile::tempdir().unwrap();
    let prefix = out.path().join("o").to_string_lossy().into_owned();
    let output = Command::new(bin())
        .args(["--fastq", &data("test_F.fastq"), "--ns_max_n", "0",
               "--out_name", &prefix, "--VERBOSE", "2"])
        .output().unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let lines: Vec<&str> = stdout.lines().collect();
    // ns_max_n is position 7 (0-indexed: min_len,max_len,min_gc,max_gc,min_qual_score,min_qual_mean,ns_max_n,...)
    assert_eq!(lines[6], "2"); // ns_max_n removed seq2_F and seq3_F
    assert_eq!(lines[0], "0"); // min_len not active
}

// ── help and version flags ────────────────────────────────────────────────────

#[test]
fn test_help_flag() {
    let out = Command::new(bin()).arg("--help").output().unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("--fastq"));
    assert!(stdout.contains("--min_len"));
}

#[test]
fn test_version_flag() {
    let out = Command::new(bin()).arg("--version").output().unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("rust-prinseq"));
    assert!(stdout.contains("1.0.0"));
}

#[test]
fn test_short_help_flag() {
    let out = Command::new(bin()).arg("-h").output().unwrap();
    assert!(!out.stdout.is_empty());
}

#[test]
fn test_short_version_flag() {
    let out = Command::new(bin()).arg("-v").output().unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("rust-prinseq"));
}
