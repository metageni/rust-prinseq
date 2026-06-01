mod reads;

use bloomfilter::Bloom;
use clap::Parser;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use rand::Rng;
use reads::{PairRead, ReadStatus, SingleRead};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    sync::{Arc, Mutex},
};

// ── CLI ───────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(name = "rust-prinseq", version, about = "PRINSEQ++ – fast FASTQ/FASTA QC filter (Rust port)")]
#[command(rename_all = "snake_case", disable_help_flag = true, disable_version_flag = true)]
struct Cli {
    /// Print help (-h)
    #[arg(short = 'h', long, action = clap::ArgAction::Help)]
    help: (),

    /// Print version (-v)
    #[arg(short = 'v', long, action = clap::ArgAction::Version)]
    version: (),

    /// Input FASTQ (or FASTA with --FASTA). Accepts .gz.
    #[arg(long)]
    fastq: String,

    /// Second input for paired-end. Accepts .gz.
    #[arg(long)]
    fastq2: Option<String>,

    /// Input is FASTA (no quality). Quality treated as 31 (A) for all bases.
    #[arg(long = "FASTA", default_value_t = false)]
    fasta: bool,

    /// Input quality is Phred+64 (old Illumina/Solexa).
    #[arg(long, default_value_t = false)]
    phred64: bool,

    /// Output format: 0 = FASTQ, 1 = FASTA.
    #[arg(long, default_value_t = 0)]
    out_format: u8,

    /// Base name for output files [default: random 6-char string].
    #[arg(long)]
    out_name: Option<String>,

    /// Write gzip-compressed output.
    #[arg(long, default_value_t = false)]
    out_gz: bool,

    /// Remove the header from the '+' line (fastq only).
    #[arg(long, default_value_t = false)]
    rm_header: bool,

    // ── per-file output overrides ──────────────────────────────────────────
    #[arg(long)] out_good:    Option<String>,
    #[arg(long)] out_good2:   Option<String>,
    #[arg(long)] out_single:  Option<String>,
    #[arg(long)] out_single2: Option<String>,
    #[arg(long)] out_bad:     Option<String>,
    #[arg(long)] out_bad2:    Option<String>,

    // ── filters ───────────────────────────────────────────────────────────
    /// Filter sequences shorter than min_len.
    #[arg(long, default_value_t = 0)]   min_len: usize,
    /// Filter sequences longer than max_len.
    #[arg(long, default_value_t = 0)]   max_len: usize,
    /// Filter sequences with GC% below min_gc.
    #[arg(long, default_value_t = 0.0)] min_gc:  f32,
    /// Filter sequences with GC% above max_gc.
    #[arg(long, default_value_t = 100.0)] max_gc: f32,
    /// Filter sequences with any base quality below min_qual_score.
    #[arg(long, default_value_t = 0)]   min_qual_score: u8,
    /// Filter sequences with mean quality below min_qual_mean.
    #[arg(long, default_value_t = 0)]   min_qual_mean:  u8,
    /// Filter sequences with more than ns_max_n Ns.
    #[arg(long, default_value_t = -1)]  ns_max_n: i32,
    /// Filter sequences with non-ACGTUN characters.
    #[arg(long, default_value_t = false)] noiupac: bool,
    /// Filter exact duplicate sequences using a bloom filter.
    #[arg(long, default_value_t = false)] derep: bool,

    /// Low-complexity entropy filter threshold [0-1] (default 0.5 if flag given without value).
    /// Scores are in range 0-1 (NOT 0-100 as in old prinseq-lite). Equivalent to ~50 in old prinseq.
    #[arg(long, num_args = 0..=1, default_missing_value = "0.5")]
    lc_entropy: Option<f32>,

    /// Low-complexity DUST filter threshold [0-1] (default 0.5 if flag given without value).
    /// Scores are in range 0-1 (NOT 0-100 as in old prinseq-lite). A score of 0.5 is roughly
    /// equivalent to a DUST score of ~20 in old prinseq-lite (i.e. use 0.2 for a lenient cutoff).
    #[arg(long, num_args = 0..=1, default_missing_value = "0.5")]
    lc_dust: Option<f32>,

    // ── trimmers ──────────────────────────────────────────────────────────
    /// Trim N bases from the 5' end.
    #[arg(long, default_value_t = 0)] trim_left:  usize,
    /// Trim N bases from the 3' end.
    #[arg(long, default_value_t = 0)] trim_right: usize,
    /// Trim poly-A/T/G tail from 5' end with minimum length N (issue #15: includes G for 2-colour Illumina).
    #[arg(long, default_value_t = 0)] trim_tail_left:  usize,
    /// Trim poly-A/T/G tail from 3' end with minimum length N (issue #15: includes G for 2-colour Illumina).
    #[arg(long, default_value_t = 0)] trim_tail_right: usize,

    /// Quality-trim 3' end; threshold (default 20 if flag given without value).
    #[arg(long, num_args = 0..=1, default_missing_value = "20")]
    trim_qual_right: Option<f32>,

    /// Quality-trim 5' end (left/5' end of read); threshold (default 20 if flag given without value).
    /// Note: 'left' refers to the 5' end of the read (issue #14).
    #[arg(long, num_args = 0..=1, default_missing_value = "20")]
    trim_qual_left: Option<f32>,

    /// Score type for quality trimming: min, mean, max, sum [default: mean].
    #[arg(long, default_value = "mean")] trim_qual_type:   String,
    /// Rule for quality trimming: lt, gt, et [default: lt].
    #[arg(long, default_value = "lt")]   trim_qual_rule:   String,
    /// Window size for quality trimming [default: 5].
    #[arg(long, default_value_t = 5)]    trim_qual_window: usize,
    /// Step size for quality trimming [default: 2].
    #[arg(long, default_value_t = 2)]    trim_qual_step:   usize,

    /// Trim sequence to at most N bases from the 5' end (issue #21).
    #[arg(long, default_value_t = 0)] trim_to_len: usize,

    /// Number of threads [default: 1].
    #[arg(long, default_value_t = 1)] threads: usize,

    /// Verbosity: 0=silent, 1=active filters only, 2=all counts one-per-line.
    #[arg(long = "VERBOSE", default_value_t = 1)]
    verbose: u8,
}

// ── helpers ───────────────────────────────────────────────────────────────────

fn random_string(n: usize) -> String {
    const CHARS: &[u8] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    let mut rng = rand::thread_rng();
    (0..n).map(|_| CHARS[rng.gen_range(0..CHARS.len())] as char).collect()
}

fn open_reader(path: &str) -> Box<dyn BufRead + Send> {
    let file = File::open(path).unwrap_or_else(|_| { eprintln!("Error: cannot open {path}"); std::process::exit(1); });
    if path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

fn open_writer(path: &str, gz: bool) -> Box<dyn Write + Send> {
    let file = File::create(path).unwrap_or_else(|_| { eprintln!("Error: cannot create {path}"); std::process::exit(1); });
    if gz {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    }
}

fn ext(out_format: u8, gz: bool) -> String {
    let base = if out_format == 1 { "fasta" } else { "fastq" };
    if gz { format!("{base}.gz") } else { base.to_string() }
}

fn write_read(read: &SingleRead, out_format: u8,
              good: &Mutex<Box<dyn Write + Send>>,
              single: &Mutex<Box<dyn Write + Send>>,
              bad: &Mutex<Box<dyn Write + Send>>) {
    let stream = match read.status {
        ReadStatus::Good   => good,
        ReadStatus::Single => single,
        ReadStatus::Bad    => bad,
    };
    let mut w = stream.lock().unwrap();
    if out_format == 1 { read.write_fasta(&mut **w); } else { read.write_fastq(&mut **w); }
}

// ── stats ─────────────────────────────────────────────────────────────────────

#[derive(Default)]
struct Stats {
    min_len: u64, max_len: u64, min_gc: u64, max_gc: u64,
    min_qual_score: u64, min_qual_mean: u64, ns_max_n: u64,
    noiupac: u64, derep: u64, lc_entropy: u64, lc_dust: u64,
    trim_tail_left: u64, trim_tail_right: u64,
    trim_qual_left: u64, trim_qual_right: u64,
    trim_left: u64, trim_right: u64, trim_to_len: u64,
}

impl Stats {
    fn add(&mut self, other: &Stats) {
        self.min_len        += other.min_len;
        self.max_len        += other.max_len;
        self.min_gc         += other.min_gc;
        self.max_gc         += other.max_gc;
        self.min_qual_score += other.min_qual_score;
        self.min_qual_mean  += other.min_qual_mean;
        self.ns_max_n       += other.ns_max_n;
        self.noiupac        += other.noiupac;
        self.derep          += other.derep;
        self.lc_entropy     += other.lc_entropy;
        self.lc_dust        += other.lc_dust;
        self.trim_tail_left += other.trim_tail_left;
        self.trim_tail_right+= other.trim_tail_right;
        self.trim_qual_left += other.trim_qual_left;
        self.trim_qual_right+= other.trim_qual_right;
        self.trim_left      += other.trim_left;
        self.trim_right     += other.trim_right;
        self.trim_to_len    += other.trim_to_len;
    }

    /// Match C++ verbose::print() exactly.
    fn print(&self, verbose: u8) {
        if verbose == 0 { return; }
        // ordered exactly as in verbose.cpp
        let fields: &[(&str, u64)] = &[
            ("min_len",         self.min_len),
            ("max_len",         self.max_len),
            ("min_gc",          self.min_gc),
            ("max_gc",          self.max_gc),
            ("min_qual_score",  self.min_qual_score),
            ("min_qual_mean",   self.min_qual_mean),
            ("ns_max_n",        self.ns_max_n),
            ("noiupac",         self.noiupac),
            ("derep",           self.derep),
            ("lc_entropy",      self.lc_entropy),
            ("lc_dust",         self.lc_dust),
            ("trim_tail_left",  self.trim_tail_left),
            ("trim_tail_right", self.trim_tail_right),
            ("trim_qual_left",  self.trim_qual_left),
            ("trim_qual_right", self.trim_qual_right),
            ("trim_left",       self.trim_left),
            ("trim_right",      self.trim_right),
        ];
        if verbose == 2 {
            // one value per line, matching C++ VERBOSE=2
            for (_, v) in fields { println!("{v}"); }
        } else {
            // VERBOSE=1: "N reads removed by -filter_name"
            for (name, v) in fields {
                if *v > 0 { println!("{v} reads removed by -{name}"); }
            }
        }
    }
}

// ── main ──────────────────────────────────────────────────────────────────────

fn main() {
    let cli = Cli::parse();

    let qual_mode: u8 = if cli.phred64 { 64 } else { 33 };
    let out_name = cli.out_name.clone().unwrap_or_else(|| random_string(6));
    let e = ext(cli.out_format, cli.out_gz);

    let is_pair = cli.fastq2.is_some();

    // ── open outputs ──────────────────────────────────────────────────────
    let mk = |name: &Option<String>, default: &str| -> Arc<Mutex<Box<dyn Write + Send>>> {
        let path = name.as_deref().unwrap_or(default);
        // compress if out_gz is set; for explicit override paths, detect by extension
        let gz = if name.is_some() { path.ends_with(".gz") } else { cli.out_gz };
        Arc::new(Mutex::new(open_writer(path, gz)))
    };

    let (good1, single1, bad1, good2, single2, bad2);
    if is_pair {
        good1   = mk(&cli.out_good,    &format!("{out_name}_good_out_R1.{e}"));
        single1 = mk(&cli.out_single,  &format!("{out_name}_single_out_R1.{e}"));
        bad1    = mk(&cli.out_bad,     &format!("{out_name}_bad_out_R1.{e}"));
        good2   = mk(&cli.out_good2,   &format!("{out_name}_good_out_R2.{e}"));
        single2 = mk(&cli.out_single2, &format!("{out_name}_single_out_R2.{e}"));
        bad2    = mk(&cli.out_bad2,    &format!("{out_name}_bad_out_R2.{e}"));
    } else {
        good1   = mk(&cli.out_good, &format!("{out_name}_good_out.{e}"));
        bad1    = mk(&cli.out_bad,  &format!("{out_name}_bad_out.{e}"));
        single1 = Arc::clone(&bad1);
        good2   = Arc::clone(&good1);
        single2 = Arc::clone(&good1);
        bad2    = Arc::clone(&bad1);
    }

    // ── bloom filter for derep ────────────────────────────────────────────
    let bloom: Option<Arc<Mutex<Bloom<String>>>> = if cli.derep {
        Some(Arc::new(Mutex::new(Bloom::new_for_fp_rate(10_000_000, 0.000001))))
    } else { None };

    // ── shared reader(s) ─────────────────────────────────────────────────
    let reader1 = Arc::new(Mutex::new(open_reader(&cli.fastq)));
    let reader2 = cli.fastq2.as_deref().map(|p| Arc::new(Mutex::new(open_reader(p))));

    let global_stats = Arc::new(Mutex::new(Stats::default()));
    let threads = cli.threads.max(1);
    let cli = Arc::new(cli);

    let handles: Vec<_> = (0..threads).map(|_| {
        let cli     = Arc::clone(&cli);
        let reader1 = Arc::clone(&reader1);
        let reader2 = reader2.as_ref().map(Arc::clone);
        let good1   = Arc::clone(&good1);
        let single1 = Arc::clone(&single1);
        let bad1    = Arc::clone(&bad1);
        let good2   = Arc::clone(&good2);
        let single2 = Arc::clone(&single2);
        let bad2    = Arc::clone(&bad2);
        let bloom   = bloom.as_ref().map(Arc::clone);
        let gstats  = Arc::clone(&global_stats);

        std::thread::spawn(move || {
            let mut stats = Stats::default();

            if is_pair {
                let r2 = reader2.unwrap();
                let mut pair = PairRead::new(qual_mode);
                loop {
                    {
                        let mut r1 = reader1.lock().unwrap();
                        let mut r2 = r2.lock().unwrap();
                        let ok1 = if cli.fasta { pair.read1.read_fasta(&mut **r1) } else { pair.read1.read_fastq(&mut **r1) };
                        let ok2 = if cli.fasta { pair.read2.read_fasta(&mut **r2) } else { pair.read2.read_fastq(&mut **r2) };
                        if !ok1 || !ok2 { break; }
                        pair.read1.status = ReadStatus::Good;
                        pair.read2.status = ReadStatus::Good;
                    }
                    apply_pair_filters(&cli, &mut pair, &bloom, &mut stats);
                    write_read(&pair.read1, cli.out_format, &good1, &single1, &bad1);
                    write_read(&pair.read2, cli.out_format, &good2, &single2, &bad2);
                }
            } else {
                let mut read = SingleRead::new(qual_mode);
                loop {
                    {
                        let mut r = reader1.lock().unwrap();
                        let ok = if cli.fasta { read.read_fasta(&mut **r) } else { read.read_fastq(&mut **r) };
                        if !ok { break; }
                        read.status = ReadStatus::Good;
                    }
                    apply_single_filters(&cli, &mut read, &bloom, &mut stats);
                    write_read(&read, cli.out_format, &good1, &single1, &bad1);
                }
            }

            gstats.lock().unwrap().add(&stats);
        })
    }).collect();

    for h in handles { h.join().unwrap(); }
    global_stats.lock().unwrap().print(cli.verbose);
}

// ── filter/trim pipelines ─────────────────────────────────────────────────────

fn apply_single_filters(
    cli: &Cli,
    read: &mut SingleRead,
    bloom: &Option<Arc<Mutex<Bloom<String>>>>,
    stats: &mut Stats,
) {
    if cli.trim_left  > 0 { stats.trim_left  += read.trim_left(cli.trim_left)   as u64; }
    if cli.trim_right > 0 { stats.trim_right += read.trim_right(cli.trim_right) as u64; }
    if cli.trim_tail_left  > 0 { stats.trim_tail_left  += read.trim_tail_left(cli.trim_tail_left)   as u64; }
    if cli.trim_tail_right > 0 { stats.trim_tail_right += read.trim_tail_right(cli.trim_tail_right) as u64; }
    if let Some(thr) = cli.trim_qual_right {
        stats.trim_qual_right += read.trim_qual_right(&cli.trim_qual_type, &cli.trim_qual_rule, cli.trim_qual_step, cli.trim_qual_window, thr) as u64;
    }
    if let Some(thr) = cli.trim_qual_left {
        stats.trim_qual_left += read.trim_qual_left(&cli.trim_qual_type, &cli.trim_qual_rule, cli.trim_qual_step, cli.trim_qual_window, thr) as u64;
    }
    if cli.trim_to_len > 0 { stats.trim_to_len += read.trim_to_len(cli.trim_to_len) as u64; }
    if cli.ns_max_n >= 0 { stats.ns_max_n += read.filter_ns_max_n(cli.ns_max_n as usize) as u64; }
    if cli.min_qual_mean  > 0 { stats.min_qual_mean  += read.filter_min_qual_mean(cli.min_qual_mean as f32)  as u64; }
    if cli.min_qual_score > 0 { stats.min_qual_score += read.filter_min_qual_score(cli.min_qual_score)       as u64; }
    if cli.noiupac        { stats.noiupac += read.filter_noiupac() as u64; }
    if cli.min_len > 0    { stats.min_len += read.filter_min_len(cli.min_len) as u64; }
    if cli.max_len > 0    { stats.max_len += read.filter_max_len(cli.max_len) as u64; }
    if cli.max_gc < 100.0 { stats.max_gc  += read.filter_max_gc(cli.max_gc)  as u64; }
    if cli.min_gc > 0.0   { stats.min_gc  += read.filter_min_gc(cli.min_gc)  as u64; }
    if let Some(b) = bloom {
        let mut bf = b.lock().unwrap();
        let dup = bf.check(&read.seq);
        if !dup { bf.set(&read.seq); }
        if dup { read.set_status(ReadStatus::Bad); stats.derep += 1; }
    }
    if let Some(thr) = cli.lc_entropy { stats.lc_entropy += read.filter_entropy(thr) as u64; }
    if let Some(thr) = cli.lc_dust    { stats.lc_dust    += read.filter_dust(thr)    as u64; }
    if cli.rm_header { read.rm_header(); }
}

fn apply_pair_filters(
    cli: &Cli,
    pair: &mut PairRead,
    bloom: &Option<Arc<Mutex<Bloom<String>>>>,
    stats: &mut Stats,
) {
    if cli.trim_left  > 0 { stats.trim_left  += pair.trim_left(cli.trim_left)   as u64; }
    if cli.trim_right > 0 { stats.trim_right += pair.trim_right(cli.trim_right) as u64; }
    if cli.trim_tail_left  > 0 { stats.trim_tail_left  += pair.trim_tail_left(cli.trim_tail_left)   as u64; }
    if cli.trim_tail_right > 0 { stats.trim_tail_right += pair.trim_tail_right(cli.trim_tail_right) as u64; }
    if let Some(thr) = cli.trim_qual_right {
        stats.trim_qual_right += pair.trim_qual_right(&cli.trim_qual_type, &cli.trim_qual_rule, cli.trim_qual_step, cli.trim_qual_window, thr) as u64;
    }
    if let Some(thr) = cli.trim_qual_left {
        stats.trim_qual_left += pair.trim_qual_left(&cli.trim_qual_type, &cli.trim_qual_rule, cli.trim_qual_step, cli.trim_qual_window, thr) as u64;
    }
    if cli.trim_to_len > 0 { stats.trim_to_len += pair.trim_to_len(cli.trim_to_len) as u64; }
    if cli.ns_max_n >= 0 { stats.ns_max_n += pair.filter_ns_max_n(cli.ns_max_n as usize) as u64; }
    if cli.min_qual_mean  > 0 { stats.min_qual_mean  += pair.filter_min_qual_mean(cli.min_qual_mean as f32)  as u64; }
    if cli.min_qual_score > 0 { stats.min_qual_score += pair.filter_min_qual_score(cli.min_qual_score)       as u64; }
    if cli.noiupac        { stats.noiupac += pair.filter_noiupac() as u64; }
    if cli.min_len > 0    { stats.min_len += pair.filter_min_len(cli.min_len) as u64; }
    if cli.max_len > 0    { stats.max_len += pair.filter_max_len(cli.max_len) as u64; }
    if cli.max_gc < 100.0 { stats.max_gc  += pair.filter_max_gc(cli.max_gc)  as u64; }
    if cli.min_gc > 0.0   { stats.min_gc  += pair.filter_min_gc(cli.min_gc)  as u64; }
    if let Some(b) = bloom {
        let mut bf = b.lock().unwrap();
        let dup1 = bf.check(&pair.read1.seq);
        let dup2 = bf.check(&pair.read2.seq);
        if !dup1 { bf.set(&pair.read1.seq); }
        if !dup2 { bf.set(&pair.read2.seq); }
        pair.set_derep_status(dup1, dup2);
        stats.derep += (dup1 as u64) + (dup2 as u64);
    }
    if let Some(thr) = cli.lc_entropy { stats.lc_entropy += pair.filter_entropy(thr) as u64; }
    if let Some(thr) = cli.lc_dust    { stats.lc_dust    += pair.filter_dust(thr)    as u64; }
    if cli.rm_header { pair.rm_header(); }
}
