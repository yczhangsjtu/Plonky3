impl<
        F,
        MDS,
        ISL,
        const WIDTH: usize,
        const CAPACITY: usize,
        const ALPHA: u64,
        const SEC_LEVEL: usize,
    > CryptographicPermutation<[F; WIDTH]>
    for Rescue<F, MDS, ISL, WIDTH, CAPACITY, ALPHA, SEC_LEVEL>
where
    F: PrimeField,
    MDS: MDSPermutation<F, WIDTH>,
    ISL: InverseSboxLayer<F, WIDTH, ALPHA>,
{
    fn permute(&self, state: [F; WIDTH]) -> [F; WIDTH] {
        // Rescue-XLIX permutation

        let mut state = state;

        for round in 0..self.num_rounds {
            // S-box
            Self::sbox_layer(&mut state);

            // MDS
            self.mds.permute_mut(&mut state);

            // Constants
            for j in 0..WIDTH {
                state[j] += self.round_constants[round * WIDTH * 2 + j];
            }

            // Inverse S-box
            self.isl.inverse_sbox_layer(&mut state);

            // MDS
            self.mds.permute_mut(&mut state);

            // Constants
            for j in 0..WIDTH {
                state[j] += self.round_constants[round * WIDTH * 2 + WIDTH + j];
            }
        }

        state
    }
}

impl<
        F,
        MDS,
        ISL,
        const WIDTH: usize,
        const CAPACITY: usize,
        const ALPHA: u64,
        const SEC_LEVEL: usize,
    > ArrayPermutation<F, WIDTH> for Rescue<F, MDS, ISL, WIDTH, CAPACITY, ALPHA, SEC_LEVEL>
where
    F: PrimeField,
    MDS: MDSPermutation<F, WIDTH>,
    ISL: InverseSboxLayer<F, WIDTH, ALPHA>,
{
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_blake3() {
        
    }
}
