use std::collections::HashMap;
use std::io::{BufRead, Write};

/// Read status: 0 = good, 1 = single (mate failed), 2 = bad
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum ReadStatus {
    Good = 0,
    Single = 1,
    Bad = 2,
}

pub struct SingleRead {
    pub name: String,
    pub seq: String,
    pub sep: String,
    pub qual: String,
    pub status: ReadStatus,
    pub qual_mode: u8, // 33 = phred33, 64 = phred64
}

impl SingleRead {
    pub fn new(qual_mode: u8) -> Self {
        SingleRead {
            name: String::new(),
            seq: String::new(),
            sep: String::new(),
            qual: String::new(),
            status: ReadStatus::Good,
            qual_mode,
        }
    }

    /// Read one record from a FASTQ or FASTA stream. Returns false at EOF.
    pub fn read_fastq(&mut self, reader: &mut dyn BufRead) -> bool {
        self.status = ReadStatus::Good;
        self.name.clear();
        self.seq.clear();
        self.sep.clear();
        self.qual.clear();

        if std::io::BufRead::read_line(reader, &mut self.name).unwrap_or(0) == 0 {
            return false;
        }
        let name = self.name.trim_end().to_string();
        if name.is_empty() {
            return false;
        }
        self.name = name;

        self.seq.clear();
        std::io::BufRead::read_line(reader, &mut self.seq).ok();
        self.seq = self.seq.trim_end().to_string();

        self.sep.clear();
        std::io::BufRead::read_line(reader, &mut self.sep).ok();
        self.sep = self.sep.trim_end().to_string();

        self.qual.clear();
        std::io::BufRead::read_line(reader, &mut self.qual).ok();
        self.qual = self.qual.trim_end().to_string();

        true
    }

    /// Read one record from a FASTA stream. Quality is set to 'A' (phred 32).
    pub fn read_fasta(&mut self, reader: &mut dyn BufRead) -> bool {
        self.status = ReadStatus::Good;
        self.name.clear();
        self.seq.clear();

        // Skip blank lines / find '>'
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line).unwrap_or(0) == 0 {
                return false;
            }
            let t = line.trim_end();
            if t.starts_with('>') {
                self.name = format!("@{}", &t[1..]);
                break;
            }
        }

        // Read sequence lines until next '>' or EOF
        loop {
            let pos = reader.fill_buf().map(|b| b.first().copied()).unwrap_or(None);
            if pos == Some(b'>') || pos.is_none() {
                break;
            }
            line.clear();
            if reader.read_line(&mut line).unwrap_or(0) == 0 {
                break;
            }
            self.seq.push_str(line.trim_end());
        }

        self.sep = "+".to_string();
        self.qual = "A".repeat(self.seq.len());
        true
    }

    pub fn set_status(&mut self, s: ReadStatus) {
        if s > self.status {
            self.status = s;
        }
    }

    // ── Filters ──────────────────────────────────────────────────────────────

    pub fn filter_ns_max_n(&mut self, max_n: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        let n = self.seq.bytes().filter(|&b| b == b'N' || b == b'n').count();
        if n > max_n { self.set_status(ReadStatus::Bad); true } else { false }
    }

    pub fn filter_min_qual_score(&mut self, min_q: u8) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        let mode = self.qual_mode;
        if self.qual.bytes().any(|b| b.saturating_sub(mode) < min_q) {
            self.set_status(ReadStatus::Bad); true
        } else { false }
    }

    pub fn filter_min_qual_mean(&mut self, min_q: f32) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        let mode = self.qual_mode as f32;
        let sum: f32 = self.qual.bytes().map(|b| b as f32 - mode).sum();
        let mean = sum / self.qual.len() as f32;
        if mean < min_q { self.set_status(ReadStatus::Bad); true } else { false }
    }

    pub fn filter_noiupac(&mut self) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if self.seq.bytes().any(|b| !matches!(b, b'A'|b'C'|b'G'|b'T'|b'U'|b'N'|b'a'|b'c'|b'g'|b't'|b'u'|b'n')) {
            self.set_status(ReadStatus::Bad); true
        } else { false }
    }

    pub fn filter_min_len(&mut self, len: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if self.seq.len() < len { self.set_status(ReadStatus::Bad); true } else { false }
    }

    pub fn filter_max_len(&mut self, len: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if self.seq.len() > len { self.set_status(ReadStatus::Bad); true } else { false }
    }

    fn gc_pct(&self) -> f32 {
        let gc = self.seq.bytes().filter(|&b| matches!(b, b'G'|b'C'|b'g'|b'c')).count();
        100.0 * gc as f32 / self.seq.len() as f32
    }

    pub fn filter_max_gc(&mut self, max_gc: f32) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if self.gc_pct() > max_gc { self.set_status(ReadStatus::Bad); true } else { false }
    }

    pub fn filter_min_gc(&mut self, min_gc: f32) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if self.gc_pct() < min_gc { self.set_status(ReadStatus::Bad); true } else { false }
    }

    /// Shannon-Wiener entropy on 64-base windows with 3-mers
    pub fn filter_entropy(&mut self, threshold: f32) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        let seq = self.seq.as_bytes();
        let mut vals: Vec<f32> = Vec::new();
        let mut j = 0usize;
        loop {
            let start = j * 32;
            if start >= seq.len() { break; }
            let window = &seq[start..std::cmp::min(start + 64, seq.len())];
            if window.len() < 15 {
                if vals.is_empty() { vals.push(0.0); }
                break;
            }
            let mut counts: HashMap<[u8; 3], u32> = HashMap::new();
            for i in 0..window.len().saturating_sub(2) {
                let tri = [window[i], window[i+1], window[i+2]];
                *counts.entry(tri).or_insert(0) += 1;
            }
            let l = (window.len() - 2) as f64;
            let k = 64f64.min(l);
            let entropy: f64 = counts.values().map(|&c| {
                let p = c as f64 / l;
                -p * (p.ln() / k.ln())
            }).sum();
            vals.push(entropy as f32);
            j += 1;
        }
        let mean = vals.iter().sum::<f32>() / vals.len() as f32;
        if mean < threshold { self.set_status(ReadStatus::Bad); true } else { false }
    }

    /// DUST score on 64-base windows with 3-mers
    pub fn filter_dust(&mut self, threshold: f32) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        let seq = self.seq.as_bytes();
        let mut vals: Vec<f32> = Vec::new();
        let mut j = 0usize;
        loop {
            let start = j * 32;
            if start >= seq.len() { break; }
            let window = &seq[start..std::cmp::min(start + 64, seq.len())];
            if window.len() < 15 {
                if vals.is_empty() { vals.push(62.0); }
                break;
            }
            let mut counts: HashMap<[u8; 3], u32> = HashMap::new();
            for i in 0..window.len().saturating_sub(2) {
                let tri = [window[i], window[i+1], window[i+2]];
                *counts.entry(tri).or_insert(0) += 1;
            }
            let l = (window.len() - 2) as f64;
            let dust: f64 = counts.values().map(|&c| {
                let cf = c as f64;
                cf * (cf - 1.0) / (l - 1.0)
            }).sum();
            vals.push(dust as f32);
            j += 1;
        }
        let mean = vals.iter().sum::<f32>() / vals.len() as f32;
        let score = (mean * 0.5) / 31.0;
        if score > threshold { self.set_status(ReadStatus::Bad); true } else { false }
    }

    // ── Trimmers ─────────────────────────────────────────────────────────────

    pub fn trim_left(&mut self, n: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if n >= self.seq.len() { self.set_status(ReadStatus::Bad); return true; }
        self.seq.drain(..n);
        self.qual.drain(..n);
        false
    }

    pub fn trim_right(&mut self, n: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if n >= self.seq.len() { self.set_status(ReadStatus::Bad); return true; }
        let new_len = self.seq.len() - n;
        self.seq.truncate(new_len);
        self.qual.truncate(new_len);
        false
    }

    pub fn trim_tail_left(&mut self, min_len: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        // Trim poly-A/T/G tails from the 5' end (issue #15: include G for 2-colour Illumina)
        let count = self.seq.bytes().take_while(|&b| matches!(b, b'A'|b'T'|b'G'|b'a'|b't'|b'g')).count();
        if count == self.seq.len() { self.set_status(ReadStatus::Bad); return true; }
        if count >= min_len {
            self.seq.drain(..count);
            self.qual.drain(..count);
        }
        false
    }

    pub fn trim_tail_right(&mut self, min_len: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        // Trim poly-A/T/G tails from the 3' end (issue #15: include G for 2-colour Illumina)
        let count = self.seq.bytes().rev().take_while(|&b| matches!(b, b'A'|b'T'|b'G'|b'a'|b't'|b'g')).count();
        if count == self.seq.len() { self.set_status(ReadStatus::Bad); return true; }
        if count >= min_len {
            let new_len = self.seq.len() - count;
            self.seq.truncate(new_len);
            self.qual.truncate(new_len);
        }
        false
    }

    /// Trim sequence to at most `len` bases from the 5' end (issue #21).
    pub fn trim_to_len(&mut self, len: usize) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        if self.seq.len() > len {
            self.seq.truncate(len);
            self.qual.truncate(len);
        }
        false
    }

    fn window_score(window: &[u8], qual_mode: u8, score_type: &str) -> f32 {
        let scores: Vec<f32> = window.iter().map(|&b| (b.saturating_sub(qual_mode)) as f32).collect();
        match score_type {
            "min"  => scores.iter().cloned().fold(f32::INFINITY, f32::min),
            "max"  => scores.iter().cloned().fold(f32::NEG_INFINITY, f32::max),
            "sum"  => scores.iter().sum(),
            _      => scores.iter().sum::<f32>() / scores.len() as f32, // mean
        }
    }

    pub fn trim_qual_right(&mut self, score_type: &str, rule: &str, step: usize, window: usize, threshold: f32) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        let mode = self.qual_mode;
        loop {
            let len = self.qual.len();
            if len < window { break; }
            let w = &self.qual.as_bytes()[len - window..];
            let score = Self::window_score(w, mode, score_type);
            let fail = match rule { "gt" => score > threshold, "et" => score == threshold, _ => score < threshold };
            if fail {
                let trim = step.min(len);
                self.seq.truncate(len - trim);
                self.qual.truncate(len - trim);
            } else { break; }
        }
        if self.qual.is_empty() { self.set_status(ReadStatus::Bad); true } else { false }
    }

    pub fn trim_qual_left(&mut self, score_type: &str, rule: &str, step: usize, window: usize, threshold: f32) -> bool {
        if self.status == ReadStatus::Bad { return false; }
        let mode = self.qual_mode;
        loop {
            let len = self.qual.len();
            if len < window { break; }
            let w = &self.qual.as_bytes()[..window];
            let score = Self::window_score(w, mode, score_type);
            let fail = match rule { "gt" => score > threshold, "et" => score == threshold, _ => score < threshold };
            if fail {
                let trim = step.min(len);
                self.seq.drain(..trim);
                self.qual.drain(..trim);
            } else { break; }
        }
        if self.qual.is_empty() { self.set_status(ReadStatus::Bad); true } else { false }
    }

    pub fn rm_header(&mut self) {
        self.sep = "+".to_string();
    }

    // ── Output ───────────────────────────────────────────────────────────────

    pub fn write_fastq(&self, out: &mut dyn Write) {
        writeln!(out, "{}\n{}\n{}\n{}", self.name, self.seq, self.sep, self.qual).ok();
    }

    pub fn write_fasta(&self, out: &mut dyn Write) {
        let header = format!(">{}", &self.name[1..]);
        writeln!(out, "{}\n{}", header, self.seq).ok();
    }
}

// ── PairRead ─────────────────────────────────────────────────────────────────

pub struct PairRead {
    pub read1: SingleRead,
    pub read2: SingleRead,
}

impl PairRead {
    pub fn new(qual_mode: u8) -> Self {
        PairRead {
            read1: SingleRead::new(qual_mode),
            read2: SingleRead::new(qual_mode),
        }
    }

    /// After applying a filter to both reads, promote the surviving mate to "single".
    fn sync_status(&mut self) {
        match (self.read1.status, self.read2.status) {
            (ReadStatus::Good, ReadStatus::Bad)   => self.read1.set_status(ReadStatus::Single),
            (ReadStatus::Bad,  ReadStatus::Good)  => self.read2.set_status(ReadStatus::Single),
            _ => {}
        }
    }

    pub fn set_derep_status(&mut self, dup1: bool, dup2: bool) {
        match (dup1, dup2) {
            (false, false) => {}
            (false, true)  => { self.read1.set_status(ReadStatus::Single); self.read2.set_status(ReadStatus::Bad); }
            (true,  false) => { self.read1.set_status(ReadStatus::Bad);    self.read2.set_status(ReadStatus::Single); }
            (true,  true)  => { self.read1.set_status(ReadStatus::Bad);    self.read2.set_status(ReadStatus::Bad); }
        }
    }

    // Delegate every filter/trim, then sync
    pub fn filter_ns_max_n(&mut self, n: usize) -> usize {
        let h = self.read1.filter_ns_max_n(n) as usize + self.read2.filter_ns_max_n(n) as usize;
        self.sync_status(); h
    }
    pub fn filter_min_qual_score(&mut self, q: u8) -> usize {
        let h = self.read1.filter_min_qual_score(q) as usize + self.read2.filter_min_qual_score(q) as usize;
        self.sync_status(); h
    }
    pub fn filter_min_qual_mean(&mut self, q: f32) -> usize {
        let h = self.read1.filter_min_qual_mean(q) as usize + self.read2.filter_min_qual_mean(q) as usize;
        self.sync_status(); h
    }
    pub fn filter_noiupac(&mut self) -> usize {
        let h = self.read1.filter_noiupac() as usize + self.read2.filter_noiupac() as usize;
        self.sync_status(); h
    }
    pub fn filter_min_len(&mut self, l: usize) -> usize {
        let h = self.read1.filter_min_len(l) as usize + self.read2.filter_min_len(l) as usize;
        self.sync_status(); h
    }
    pub fn filter_max_len(&mut self, l: usize) -> usize {
        let h = self.read1.filter_max_len(l) as usize + self.read2.filter_max_len(l) as usize;
        self.sync_status(); h
    }
    pub fn filter_max_gc(&mut self, v: f32) -> usize {
        let h = self.read1.filter_max_gc(v) as usize + self.read2.filter_max_gc(v) as usize;
        self.sync_status(); h
    }
    pub fn filter_min_gc(&mut self, v: f32) -> usize {
        let h = self.read1.filter_min_gc(v) as usize + self.read2.filter_min_gc(v) as usize;
        self.sync_status(); h
    }
    pub fn filter_entropy(&mut self, t: f32) -> usize {
        let h = self.read1.filter_entropy(t) as usize + self.read2.filter_entropy(t) as usize;
        self.sync_status(); h
    }
    pub fn filter_dust(&mut self, t: f32) -> usize {
        let h = self.read1.filter_dust(t) as usize + self.read2.filter_dust(t) as usize;
        self.sync_status(); h
    }
    pub fn trim_left(&mut self, n: usize) -> usize {
        let h = self.read1.trim_left(n) as usize + self.read2.trim_left(n) as usize;
        self.sync_status(); h
    }
    pub fn trim_right(&mut self, n: usize) -> usize {
        let h = self.read1.trim_right(n) as usize + self.read2.trim_right(n) as usize;
        self.sync_status(); h
    }
    pub fn trim_tail_left(&mut self, n: usize) -> usize {
        let h = self.read1.trim_tail_left(n) as usize + self.read2.trim_tail_left(n) as usize;
        self.sync_status(); h
    }
    pub fn trim_tail_right(&mut self, n: usize) -> usize {
        let h = self.read1.trim_tail_right(n) as usize + self.read2.trim_tail_right(n) as usize;
        self.sync_status(); h
    }
    pub fn trim_qual_right(&mut self, st: &str, rule: &str, step: usize, win: usize, thr: f32) -> usize {
        let h = self.read1.trim_qual_right(st, rule, step, win, thr) as usize
              + self.read2.trim_qual_right(st, rule, step, win, thr) as usize;
        self.sync_status(); h
    }
    pub fn trim_qual_left(&mut self, st: &str, rule: &str, step: usize, win: usize, thr: f32) -> usize {
        let h = self.read1.trim_qual_left(st, rule, step, win, thr) as usize
              + self.read2.trim_qual_left(st, rule, step, win, thr) as usize;
        self.sync_status(); h
    }
    pub fn trim_to_len(&mut self, len: usize) -> usize {
        let h = self.read1.trim_to_len(len) as usize + self.read2.trim_to_len(len) as usize;
        self.sync_status(); h
    }
    pub fn rm_header(&mut self) {
        self.read1.rm_header();
        self.read2.rm_header();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_read(seq: &str, qual: &str) -> SingleRead {
        let mut r = SingleRead::new(33);
        r.name = "@test".into();
        r.seq  = seq.into();
        r.sep  = "+".into();
        r.qual = qual.into();
        r
    }

    fn fastq(name: &str, seq: &str, qual: &str) -> String {
        format!("{name}\n{seq}\n+\n{qual}\n")
    }

    // ── read_fastq ────────────────────────────────────────────────────────────

    #[test]
    fn read_fastq_parses_record() {
        let input = fastq("@seq1", "ACGT", "IIII");
        let mut r = SingleRead::new(33);
        assert!(r.read_fastq(&mut Cursor::new(input)));
        assert_eq!(r.name, "@seq1");
        assert_eq!(r.seq,  "ACGT");
        assert_eq!(r.qual, "IIII");
    }

    #[test]
    fn read_fastq_returns_false_at_eof() {
        let mut r = SingleRead::new(33);
        assert!(!r.read_fastq(&mut Cursor::new("")));
    }

    #[test]
    fn read_fasta_parses_record() {
        let input = ">seq1\nACGT\n";
        let mut r = SingleRead::new(33);
        assert!(r.read_fasta(&mut Cursor::new(input)));
        assert_eq!(r.name, "@seq1");
        assert_eq!(r.seq,  "ACGT");
        assert_eq!(r.qual, "AAAA"); // 'A' repeated
    }

    // ── filters ───────────────────────────────────────────────────────────────

    #[test]
    fn filter_min_len_pass() {
        let mut r = make_read("ACGTACGT", "IIIIIIII");
        assert!(!r.filter_min_len(8));
        assert_eq!(r.status, ReadStatus::Good);
    }

    #[test]
    fn filter_min_len_fail() {
        let mut r = make_read("ACGT", "IIII");
        assert!(r.filter_min_len(5));
        assert_eq!(r.status, ReadStatus::Bad);
    }

    #[test]
    fn filter_max_len_fail() {
        let mut r = make_read("ACGTACGT", "IIIIIIII");
        assert!(r.filter_max_len(4));
        assert_eq!(r.status, ReadStatus::Bad);
    }

    #[test]
    fn filter_ns_max_n_pass() {
        let mut r = make_read("ACGN", "IIII");
        assert!(!r.filter_ns_max_n(1));
    }

    #[test]
    fn filter_ns_max_n_fail() {
        let mut r = make_read("ACNN", "IIII");
        assert!(r.filter_ns_max_n(1));
        assert_eq!(r.status, ReadStatus::Bad);
    }

    #[test]
    fn filter_min_qual_score_fail() {
        // qual 'A' = ASCII 65, phred33 = 32; min=33 → fail
        let mut r = make_read("ACGT", "AAAA");
        assert!(r.filter_min_qual_score(33));
        assert_eq!(r.status, ReadStatus::Bad);
    }

    #[test]
    fn filter_min_qual_score_pass() {
        // qual 'I' = ASCII 73, phred33 = 40 ≥ 33
        let mut r = make_read("ACGT", "IIII");
        assert!(!r.filter_min_qual_score(33));
    }

    #[test]
    fn filter_min_qual_mean_fail() {
        // qual '!' = ASCII 33, phred33 = 0; mean=0 < 20 → fail
        let mut r = make_read("ACGT", "!!!!");
        assert!(r.filter_min_qual_mean(20.0));
    }

    #[test]
    fn filter_noiupac_pass() {
        let mut r = make_read("ACGTUN", "IIIIII");
        assert!(!r.filter_noiupac());
    }

    #[test]
    fn filter_noiupac_fail() {
        let mut r = make_read("ACGTX", "IIIII");
        assert!(r.filter_noiupac());
        assert_eq!(r.status, ReadStatus::Bad);
    }

    #[test]
    fn filter_gc_bounds() {
        // GGCC = 100% GC
        let mut r = make_read("GGCC", "IIII");
        assert!(r.filter_max_gc(50.0));

        let mut r2 = make_read("AATT", "IIII");
        assert!(r2.filter_min_gc(10.0));
    }

    // ── trimmers ─────────────────────────────────────────────────────────────

    #[test]
    fn trim_left_basic() {
        let mut r = make_read("ACGTACGT", "12345678");
        assert!(!r.trim_left(3));
        assert_eq!(r.seq,  "TACGT");
        assert_eq!(r.qual, "45678");
    }

    #[test]
    fn trim_left_entire_seq_becomes_bad() {
        let mut r = make_read("ACG", "III");
        assert!(r.trim_left(3));
        assert_eq!(r.status, ReadStatus::Bad);
    }

    #[test]
    fn trim_right_basic() {
        let mut r = make_read("ACGTACGT", "12345678");
        assert!(!r.trim_right(3));
        assert_eq!(r.seq,  "ACGTA");
        assert_eq!(r.qual, "12345");
    }

    #[test]
    fn trim_tail_left_removes_at() {
        // AAACGT → 3 A's at start ≥ min_len=3 → trim
        let mut r = make_read("AAACGT", "IIIIII");
        assert!(!r.trim_tail_left(3));
        assert_eq!(r.seq, "CGT");
    }

    #[test]
    fn trim_tail_left_below_min_len_no_trim() {
        // AA at start < min_len=3 → no trim
        let mut r = make_read("AACGT", "IIIII");
        assert!(!r.trim_tail_left(3));
        assert_eq!(r.seq, "AACGT");
    }

    #[test]
    fn trim_tail_right_removes_at() {
        // CGTTTT: 4 T's at end ≥ min_len=3 → trimmed. Note: G is also a tail char,
        // so "GTTTT" (5 chars) are trimmed → "C"
        let mut r = make_read("CGTTTT", "IIIIII");
        assert!(!r.trim_tail_right(3));
        assert_eq!(r.seq, "C");
    }

    // issue #15: poly-G trimming for 2-colour Illumina
    #[test]
    fn trim_tail_right_removes_poly_g() {
        // ACGGGG: 4 G's at end ≥ min_len=3 → trimmed → "AC"
        let mut r = make_read("ACGGGG", "IIIIII");
        assert!(!r.trim_tail_right(3));
        assert_eq!(r.seq, "AC");
    }

    #[test]
    fn trim_tail_left_removes_poly_g() {
        // GGGCGT: 3 G's at start ≥ min_len=3 → trimmed → "CGT"
        let mut r = make_read("GGGCGT", "IIIIII");
        assert!(!r.trim_tail_left(3));
        assert_eq!(r.seq, "CGT");
    }

    // issue #21: trim_to_len
    #[test]
    fn trim_to_len_truncates() {
        let mut r = make_read("ACGTACGT", "IIIIIIII");
        assert!(!r.trim_to_len(4));
        assert_eq!(r.seq,  "ACGT");
        assert_eq!(r.qual, "IIII");
    }

    #[test]
    fn trim_to_len_no_op_when_shorter() {
        let mut r = make_read("ACGT", "IIII");
        assert!(!r.trim_to_len(10));
        assert_eq!(r.seq, "ACGT");
    }

    #[test]
    fn trim_qual_right_trims_low_quality() {
        // qual: IIIAAA → last window of 3 = AAA, phred 32 < 40 → trim
        let mut r = make_read("ACGTAC", "IIIAAA");
        r.trim_qual_right("mean", "lt", 3, 3, 40.0);
        assert_eq!(r.seq,  "ACG");
        assert_eq!(r.qual, "III");
    }

    #[test]
    fn trim_qual_left_trims_low_quality() {
        // qual: AAAIII → first window of 3 = AAA, phred 32 < 40 → trim
        let mut r = make_read("ACGTAC", "AAAIII");
        r.trim_qual_left("mean", "lt", 3, 3, 40.0);
        assert_eq!(r.seq,  "TAC");
        assert_eq!(r.qual, "III");
    }

    #[test]
    fn trim_qual_right_empty_becomes_bad() {
        let mut r = make_read("ACG", "AAA");
        r.trim_qual_right("mean", "lt", 3, 3, 40.0);
        assert_eq!(r.status, ReadStatus::Bad);
    }

    #[test]
    fn rm_header_clears_sep() {
        let mut r = make_read("ACGT", "IIII");
        r.sep = "+seq1".into();
        r.rm_header();
        assert_eq!(r.sep, "+");
    }

    // ── set_status monotonicity ───────────────────────────────────────────────

    #[test]
    fn set_status_only_worsens() {
        let mut r = make_read("A", "I");
        r.set_status(ReadStatus::Single);
        assert_eq!(r.status, ReadStatus::Single);
        r.set_status(ReadStatus::Good); // should not improve
        assert_eq!(r.status, ReadStatus::Single);
        r.set_status(ReadStatus::Bad);
        assert_eq!(r.status, ReadStatus::Bad);
    }

    // ── skip filters when already bad ────────────────────────────────────────

    #[test]
    fn filters_skip_when_bad() {
        let mut r = make_read("ACGT", "IIII");
        r.status = ReadStatus::Bad;
        assert!(!r.filter_min_len(100)); // would fail but returns false (already bad)
        assert!(!r.filter_noiupac());
    }

    // ── PairRead sync_status ──────────────────────────────────────────────────

    #[test]
    fn pair_sync_status_promotes_to_single() {
        let mut p = PairRead::new(33);
        p.read1.seq = "ACGT".into(); p.read1.qual = "IIII".into(); p.read1.name = "@r1".into(); p.read1.sep = "+".into();
        p.read2.seq = "ACGT".into(); p.read2.qual = "IIII".into(); p.read2.name = "@r2".into(); p.read2.sep = "+".into();
        // Manually fail read2
        p.read2.status = ReadStatus::Bad;
        p.filter_min_len(1); // triggers sync
        assert_eq!(p.read1.status, ReadStatus::Single);
        assert_eq!(p.read2.status, ReadStatus::Bad);
    }
}
