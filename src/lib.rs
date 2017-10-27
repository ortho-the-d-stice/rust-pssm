
extern crate ndarray;
use std::convert::From;
use ndarray::prelude::{Array, Array2, Axis};
use std::f32;

struct Base([f32; 4]);

impl Base {
    const A: [f32; 4] = [1.0, 0.0, 0.0, 0.0];
    const T: [f32; 4] = [0.0, 1.0, 0.0, 0.0];
    const G: [f32; 4] = [0.0, 0.0, 1.0, 0.0];
    const C: [f32; 4] = [0.0, 0.0, 0.0, 1.0];

    pub fn new(seq: &[u8]) -> Array2<f32> {
        let mut flat: Vec<f32> = Vec::with_capacity(seq.len() * 4);
        for b in seq.iter() {
            flat.extend(
                match *b {
                    b'A' => Base::A,
                    b'T' => Base::T,
                    b'G' => Base::G,
                    b'C' => Base::C,
                    _ => panic!("bad base: {}", b),
                }.iter(),
            );
        }
        Array::from_shape_vec((seq.len(), 4), flat).expect("ndarray init")
    }
}

pub enum BasePos {
    A = 0,
    T = 1,
    G = 2,
    C = 3,
}
impl BasePos {
    pub fn get(base: u8) -> usize {
        match base {
            b'A' => BasePos::A as usize,
            b'T' => BasePos::T as usize,
            b'G' => BasePos::G as usize,
            b'C' => BasePos::C as usize,
            _ => panic!("unrecognized base"),
        }
    }
    pub fn put(code: usize) -> u8 {
        match code {
            0 => b'A',
            1 => b'T',
            2 => b'G',
            3 => b'C',
            _ => panic!("unrecognized code"),
        }
    }
}

/// represents motif score
#[derive(Debug)]
pub struct ScoredPos {
    pub loc: usize,
    pub sum: f32,
    pub scores: Vec<f32>,
}

impl ScoredPos {
    pub fn nil() -> ScoredPos {
        ScoredPos {
            loc: 0,
            sum: f32::NEG_INFINITY,
            scores: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Motif {
    pub seqs: Option<Vec<Vec<u8>>>,
    pub counts: Option<Array2<u16>>,
    pub pseudoct: [f32; 4],
    pub seq_ct: usize,
    pub scores: Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score: f32,
    /// sum of "best" base at each position
    pub max_score: f32,
}

impl Motif {
    /// set pseudocount for all bases
    pub fn pseudocts(&mut self, ct: f32) -> &mut Self {
        self.pseudoct = [ct; 4];
        self.calc_scores();
        self
    }

    /// set pseudocount for one base
    pub fn pseudoct(&mut self, base: u8, ct: f32) -> &mut Self {
        self.pseudoct[BasePos::get(base)] = ct;
        self.calc_scores();
        self
    }

    /// update scores field based on counts & pseudocounts
    pub fn calc_scores(&mut self) {
        match self.seqs {
            Some(_) => {
                let seqs = self.seqs.as_ref().expect("calc_scores called without seqs");
                for seq_i in 0..seqs[0].len() {
                    let mut tot: f32 = 0.0;
                    // FIXME: slices should be cleaner
                    for base_i in 0..4 {
                        tot += self.counts.as_ref().expect(
                            "calc_scores called without counts",
                        )
                            [[seq_i, base_i]] as f32;
                        tot += self.pseudoct[base_i];
                    }
                    let mut scores = self.scores.subview_mut(Axis(0), seq_i);
                    for base_i in 0..4 {
                        scores[base_i] =
                            (self.counts.as_ref().expect(
                                "calc_scores called without counts",
                            )
                                 [[seq_i, base_i]] as f32 +
                                 self.pseudoct[base_i]) / tot;
                    }
                }
            }
            _ => (),
        }
        self.calc_minmax();
    }

    pub fn calc_minmax(&mut self) {
        let pwm_len = self.len();

        // score corresponding to sum of "worst" bases at each position
        // FIXME: iter ...
        self.min_score = 0.0;
        for i in 0..pwm_len {
            // can't use the regular min/max on f32, so we use f32::min
            let min_sc = (0..4).map(|b| self.scores[[i, b]]).fold(f32::INFINITY, f32::min);
            self.min_score += min_sc;
        }

        // score corresponding to "best" base at each position
        self.max_score = 0.0;
        for i in 0..pwm_len {
            let max_sc = (0..4).map(|b| self.scores[[i, b]]).fold(f32::NEG_INFINITY, f32::max);
            self.max_score += max_sc;
        }
    }

    /// calculate consensus sequence
    pub fn consensus(&self) -> Vec<u8> {
        Vec::new()
    }

    pub fn degenerate_consensus(&self) -> Vec<u8> {
        Vec::new()
    }

    /// apply PSM to sequence, finding the offset with the highest score
    /// return None if sequence is too short
    /// see:
    ///   MATCHTM: a tool for searching transcription factor binding sites in DNA sequences
    ///   Nucleic Acids Res. 2003 Jul 1; 31(13): 3576–3579
    ///   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169193/
    ///
    pub fn score(&self, seq: &[u8]) -> Option<ScoredPos> {
        let pwm_len = self.len();
        if seq.len() < pwm_len {
            return None;
        }
        if self.max_score == self.min_score {
            return None;
        }

        let mut best_start = 0;
        let mut best_score = -1.0;
        let mut best_m = Vec::new();
        for start in 0..seq.len() - pwm_len + 1 {
            let m: Vec<f32> = (0..pwm_len).map(|i| self.scores[[i, BasePos::get(seq[start + i])]]).collect();
            let tot = m.iter().sum();
            if tot > best_score {
                best_score = tot;
                best_start = start;
                best_m = m;
            }
        }
        if best_score < 0.0 || best_score - self.min_score < 0.0 || self.max_score - self.min_score < 0.0 {
            println!("@@ found problem: best = {}; num = {} - {} = {}; demon = {}",
                     best_score, best_score, self.min_score, best_score - self.min_score, self.max_score - self.min_score);
        }
        Some(ScoredPos {
            loc: best_start,
            sum: (best_score - self.min_score) / (self.max_score - self.min_score),
            scores: best_m,
        })
    }

    pub fn len(&self) -> usize {
        self.scores.dim().0
    }
}

impl From<Vec<Vec<u8>>> for Motif {
    fn from(seqs: Vec<Vec<u8>>) -> Self {
        // null case
        if seqs.len() == 0 {
            return Motif {
                seq_ct: 0,
                seqs: None,
                counts: None,
                pseudoct: [0.0; 4],
                scores: Array2::zeros((0, 0)),
                min_score: 0.0,
                max_score: 0.0,
            };
        }

        let seqlen = seqs[0].len();
        let mut counts = Array2::zeros((seqlen, 4));

        for seq in seqs.iter() {
            if seq.len() != seqlen {
                panic!("inconsistent sequence lengths");
            }

            for (idx, base) in seq.iter().enumerate() {
                counts[[idx, BasePos::get(*base)]] += 1;
            }
        }

        let mut m = Motif {
            seq_ct: seqs.len(),
            seqs: Some(seqs),
            counts: Some(counts),
            pseudoct: [0.5; 4],
            scores: Array2::zeros((seqlen, 4)),
            min_score: 0.0,
            max_score: 0.0,
        };
        m.calc_scores();
        m
    }
}

impl From<Array2<f32>> for Motif {
    fn from(scores: Array2<f32>) -> Self {
        let mut m = Motif {
            seq_ct: 0,
            seqs: None,
            counts: None,
            pseudoct: [0.0; 4],
            scores: scores,
            min_score: 0.0,
            max_score: 0.0,
        };
        m.calc_minmax();
        m
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn simple_pwm() {
        let pwm = Motif::from(vec![
            b"AAAA".to_vec(),
            b"TTTT".to_vec(),
            b"GGGG".to_vec(),
            b"CCCC".to_vec(),
        ]);
        assert_eq!(pwm.scores, Array2::from_elem((4, 4), 0.25));
    }
    #[test]
    fn find_motif() {
        let pwm = Motif::from(vec![b"ATGC".to_vec()]);
        let seq = b"GGGGATGCGGGG";
        if let Some(ScoredPos { ref loc, ref sum, .. }) = pwm.score(seq) {
            assert_eq!(*loc, 4);
            assert_eq!(*sum, 1.0);
        } else {
            assert!(false);
        }
    }
}
