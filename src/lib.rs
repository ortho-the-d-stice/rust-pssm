extern crate ndarray;

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
}

pub struct Motif {
    pub seqs:      Vec<Vec<u8>>,
    pub counts:    Array2<u16>,
    pub pseudoct:  [f32; 4],
    seq_ct:        usize,
    scores:        Array2<f32>,
}

impl Motif {
    pub fn new(seqs: Vec<Vec<u8>>) -> Motif {
        // null case
        if seqs.len() == 0 {
            return Motif {
                seq_ct:   0,
                seqs:     Vec::new(),
                counts:   Array2::zeros((0,0)),
                pseudoct: [0.0; 4],
                scores:   Array2::zeros((0,0)),
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
            seq_ct:   seqs.len(),
            seqs:     seqs,
            counts:   counts,
            pseudoct: [0.5; 4],
            scores:   Array2::zeros((seqlen,4)),
        };
        m.calc_scores();
        m
    }

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
    fn calc_scores(&mut self) {
        for seq_i in 0 .. self.seqs[0].len() {
            let mut scores = self.scores.subview_mut(Axis(0), seq_i);
            for base_i in 0 .. 4 {
                scores[base_i] = (self.counts[[seq_i, base_i]] as f32 + self.pseudoct[base_i])
                    / self.seq_ct as f32;
            }
        }
    }

    /// calculate consensus sequence
    pub fn consensus(&self) -> Vec<u8> {
        Vec::new()
    }

    pub fn degenerate_consensus(&self) -> Vec<u8> {
        Vec::new()
    }

}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
    }
}
