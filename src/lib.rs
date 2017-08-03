
extern crate ndarray;
use std::convert::From;
use ndarray::prelude::{ Array2, Axis };

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

#[derive(Debug, Clone)]
pub struct Motif {
    pub seqs:      Option<Vec<Vec<u8>>>,
    pub counts:    Option<Array2<u16>>,
    pub pseudoct:  [f32; 4],
    pub seq_ct:        usize,
    pub scores:        Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score:     f32,
    /// sum of "best" base at each position
    pub max_score:     f32,
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
        self.pseudoct[ BasePos::get(base) ] = ct;
        self.calc_scores();
        self
    }

    /// update scores field based on counts & pseudocounts
    pub fn calc_scores(&mut self) {
        match self.seqs {
            Some(_) => {
                let seqs = self.seqs.as_ref().expect("calc_scores called without seqs");
                for seq_i in 0 .. seqs[0].len() {
                    let mut tot: f32 = 0.0;
                    // FIXME: slices should be cleaner
                    for base_i in 0 .. 4 {
                        tot += self.counts.as_ref()
                            .expect("calc_scores called without counts")[[seq_i, base_i]]
                            as f32;
                        tot += self.pseudoct[base_i];
                    }
                    let mut scores = self.scores.subview_mut(Axis(0), seq_i);
                    for base_i in 0 .. 4 {
                        scores[base_i] = (self.counts.as_ref()
                                          .expect("calc_scores called without counts")
                                          [[seq_i, base_i]] as f32 + self.pseudoct[base_i])
                            / tot;
                    }
                }
            },
            _ => ()
        }

        let (pwm_len, _) = self.scores.dim();
        // score corresponding to "worst" base at each position
        // FIXME: iter ...
        for i in 0 .. pwm_len {
            let mut min_sc = 999.9;
            for b in 0 .. 4 {
                if self.scores[[i,b]] < min_sc {
                    min_sc = self.scores[[i,b]];
                }
            }
            self.min_score += min_sc;
        }

        // score corresponding to "best" base at each position
        for i in 0 .. pwm_len {
            let mut max_sc = -999.9;
            for b in 0 .. 4 {
                if self.scores[[i,b]] > max_sc {
                    max_sc = self.scores[[i,b]];
                }
            }
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
    pub fn score(&self, seq: &[u8]) -> Option<(usize, f32)> {
        let (pwm_len, _) = self.scores.dim();
        if seq.len() < pwm_len {
            return None
        }

        let mut best_start = 0;
        let mut best_score = -1.0;
        for start in 0 .. seq.len() - pwm_len {
            let mut tot = 0.0;
            for i in 0 .. pwm_len {
                tot += self.scores[[i, BasePos::get(seq[start+i])]];
            }
            if tot > best_score {
                best_score = tot;
                best_start = start;
            }
        }

        Some((best_start, (best_score - self.min_score) / (self.max_score - self.min_score)))
    }
}

impl From<Vec<Vec<u8>>> for Motif {
    fn from(seqs: Vec<Vec<u8>>) -> Self {
        // null case
        if seqs.len() == 0 {
            return Motif {
                seq_ct:    0,
                seqs:      None,
                counts:    None,
                pseudoct:  [0.0; 4],
                scores:    Array2::zeros((0,0)),
                min_score: 0.0,
                max_score: 0.0,
            }
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
            seq_ct:    seqs.len(),
            seqs:      Some(seqs),
            counts:    Some(counts),
            pseudoct:  [0.5; 4],
            scores:    Array2::zeros((seqlen,4)),
            min_score: 0.0,
            max_score: 0.0,
        };
        m.calc_scores();
        m
    }
}

impl From<Array2<f32>> for Motif {
    fn from(scores: Array2<f32>) -> Self {
        Motif {
            seq_ct:    0,
            seqs:      None,
            counts:    None,
            pseudoct:  [0.0; 4],
            scores:    scores,
            min_score: 0.0,
            max_score: 0.0,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn simple_pwm() {
        let pwm = Motif::from( vec![ b"AAAA".to_vec(),
                                    b"TTTT".to_vec(),
                                    b"GGGG".to_vec(),
                                    b"CCCC".to_vec() ] );
        assert_eq!( pwm.scores, Array2::from_elem((4,4), 0.25) );
    }
    #[test]
    fn find_motif() {
        let pwm = Motif::from( vec![b"ATGC".to_vec()] );
        let seq = b"GGGGATGCGGGG";
        if let Some((start,score)) = pwm.score(seq) {
            assert_eq!(start, 4);
            assert_eq!(score, 1.0);
        } else {
            assert!(false);
        }
    }
}
