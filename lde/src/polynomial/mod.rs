pub(crate) mod division;

use p3_field::PrimeField32;
use p3_field::TwoAdicField;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::max;
use core::iter::Sum;
use core::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use anyhow::{ensure, Result};
use itertools::Itertools;
use p3_util::log2_strict_usize;

use crate::fft::{fft, fft_with_options, ifft, FftRootTable};
use p3_field::Field;

/// A polynomial in point-value form.
///
/// The points are implicitly `g^i`, where `g` generates the subgroup whose size equals the number
/// of points.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PolynomialValues<F: Field + PrimeField32 + TwoAdicField> {
    pub values: Vec<F>,
}

impl<F: Field + TwoAdicField + PrimeField32> PolynomialValues<F> {
    pub fn new(values: Vec<F>) -> Self {
        // Check that a subgroup exists of this size, which should be a power of two.
        debug_assert!(log2_strict_usize(values.len()) <= F::TWO_ADICITY);
        PolynomialValues { values }
    }

    pub fn constant(value: F, len: usize) -> Self {
        Self::new(vec![value; len])
    }

    pub fn zero(len: usize) -> Self {
        Self::constant(F::ZERO, len)
    }

    pub fn is_zero(&self) -> bool {
        self.values.iter().all(|x| x.is_zero())
    }

    /// Returns the polynomial whole value is one at the given index, and zero elsewhere.
    pub fn selector(len: usize, index: usize) -> Self {
        let mut result = Self::zero(len);
        result.values[index] = F::ONE;
        result
    }

    /// The number of values stored.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn ifft(self) -> PolynomialCoeffs<F> {
        ifft(self)
    }

    /// Returns the polynomial whose evaluation on the coset `shift*H` is `self`.
    pub fn coset_ifft(self, shift: F) -> PolynomialCoeffs<F> {
        let mut shifted_coeffs = self.ifft();
        shifted_coeffs
            .coeffs
            .iter_mut()
            .zip(shift.inverse().powers())
            .for_each(|(c, r)| {
                *c *= r;
            });
        shifted_coeffs
    }

    pub fn lde_multiple(polys: Vec<Self>, rate_bits: usize) -> Vec<Self> {
        polys.into_iter().map(|p| p.lde(rate_bits)).collect()
    }

    pub fn lde(self, rate_bits: usize) -> Self {
        let coeffs = ifft(self).lde(rate_bits);
        fft_with_options(coeffs, Some(rate_bits), None)
    }

    /// Low-degree extend `Self` (seen as evaluations over the subgroup) onto a coset.
    pub fn lde_onto_coset(self, rate_bits: usize) -> Self {
        let coeffs = ifft(self).lde(rate_bits);
        coeffs.coset_fft_with_options(F::coset_shift(), Some(rate_bits), None)
    }

    pub fn degree(&self) -> usize {
        self.degree_plus_one()
            .checked_sub(1)
            .expect("deg(0) is undefined")
    }

    pub fn degree_plus_one(&self) -> usize {
        self.clone().ifft().degree_plus_one()
    }

    /// Adds `rhs * rhs_weight` to `self`. Assumes `self.len() == rhs.len()`.
    pub fn add_assign_scaled(&mut self, rhs: &Self, rhs_weight: F) {
        self.values
            .iter_mut()
            .zip_eq(&rhs.values)
            .for_each(|(self_v, rhs_v)| *self_v += *rhs_v * rhs_weight)
    }
}

impl<F: Field + TwoAdicField + p3_field::PrimeField32> From<Vec<F>> for PolynomialValues<F> {
    fn from(values: Vec<F>) -> Self {
        Self::new(values)
    }
}

/// A polynomial in coefficient form.
#[derive(Clone, Debug)]
pub struct PolynomialCoeffs<F: Field + TwoAdicField + PrimeField32> {
    pub coeffs: Vec<F>,
}

impl<F: Field + TwoAdicField + PrimeField32> PolynomialCoeffs<F> {
    pub fn new(coeffs: Vec<F>) -> Self {
        PolynomialCoeffs { coeffs }
    }

    /// The empty list of coefficients, which is the smallest encoding of the zero polynomial.
    pub fn empty() -> Self {
        Self::new(Vec::new())
    }

    pub fn zero(len: usize) -> Self {
        Self::new(vec![F::ZERO; len])
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|x| x.is_zero())
    }

    /// The number of coefficients. This does not filter out any zero coefficients, so it is not
    /// necessarily related to the degree.
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn log_len(&self) -> usize {
        log2_strict_usize(self.len())
    }

    pub fn chunks(&self, chunk_size: usize) -> Vec<Self> {
        self.coeffs
            .chunks(chunk_size)
            .map(|chunk| PolynomialCoeffs::new(chunk.to_vec()))
            .collect()
    }

    pub fn eval(&self, x: F) -> F {
        self.coeffs
            .iter()
            .rev()
            .fold(F::ZERO, |acc, &c| acc * x + c)
    }

    /// Evaluate the polynomial at a point given its powers. The first power is the point itself, not 1.
    pub fn eval_with_powers(&self, powers: &[F]) -> F {
        debug_assert_eq!(self.coeffs.len(), powers.len() + 1);
        let acc = self.coeffs[0];
        self.coeffs[1..]
            .iter()
            .zip(powers)
            .fold(acc, |acc, (&x, &c)| acc + c * x)
    }

    pub fn lde_multiple(polys: Vec<&Self>, rate_bits: usize) -> Vec<Self> {
        polys.into_iter().map(|p| p.lde(rate_bits)).collect()
    }

    pub fn lde(&self, rate_bits: usize) -> Self {
        self.padded(self.len() << rate_bits)
    }

    pub fn pad(&mut self, new_len: usize) -> Result<()> {
        ensure!(
            new_len >= self.len(),
            "Trying to pad a polynomial of length {} to a length of {}.",
            self.len(),
            new_len
        );
        self.coeffs.resize(new_len, F::ZERO);
        Ok(())
    }

    pub fn padded(&self, new_len: usize) -> Self {
        let mut poly = self.clone();
        poly.pad(new_len).unwrap();
        poly
    }

    /// Removes any leading zero coefficients.
    pub fn trim(&mut self) {
        self.coeffs.truncate(self.degree_plus_one());
    }

    /// Removes some leading zero coefficients, such that a desired length is reached. Fails if a
    /// nonzero coefficient is encountered before then.
    pub fn trim_to_len(&mut self, len: usize) -> Result<()> {
        ensure!(self.len() >= len);
        ensure!(self.coeffs[len..].iter().all(F::is_zero));
        self.coeffs.truncate(len);
        Ok(())
    }

    /// Removes any leading zero coefficients.
    pub fn trimmed(&self) -> Self {
        let coeffs = self.coeffs[..self.degree_plus_one()].to_vec();
        Self { coeffs }
    }

    /// Degree of the polynomial + 1, or 0 for a polynomial with no non-zero coefficients.
    pub fn degree_plus_one(&self) -> usize {
        (0usize..self.len())
            .rev()
            .find(|&i| !self.coeffs[i].is_zero())
            .map_or(0, |i| i + 1)
    }

    /// Leading coefficient.
    pub fn lead(&self) -> F {
        self.coeffs
            .iter()
            .rev()
            .find(|x| !x.is_zero())
            .map_or(F::ZERO, |x| *x)
    }

    /// Reverse the order of the coefficients, not taking into account the leading zero coefficients.
    pub(crate) fn rev(&self) -> Self {
        Self::new(self.trimmed().coeffs.into_iter().rev().collect())
    }

    pub fn fft(self) -> PolynomialValues<F> {
        fft(self)
    }

    pub fn fft_with_options(
        self,
        zero_factor: Option<usize>,
        root_table: Option<&FftRootTable<F>>,
    ) -> PolynomialValues<F> {
        fft_with_options(self, zero_factor, root_table)
    }

    /// Returns the evaluation of the polynomial on the coset `shift*H`.
    pub fn coset_fft(&self, shift: F) -> PolynomialValues<F> {
        self.coset_fft_with_options(shift, None, None)
    }

    /// Returns the evaluation of the polynomial on the coset `shift*H`.
    pub fn coset_fft_with_options(
        &self,
        shift: F,
        zero_factor: Option<usize>,
        root_table: Option<&FftRootTable<F>>,
    ) -> PolynomialValues<F> {
        let modified_poly: Self = shift
            .powers()
            .zip(&self.coeffs)
            .map(|(r, &c)| r * c)
            .collect::<Vec<_>>()
            .into();
        modified_poly.fft_with_options(zero_factor, root_table)
    }
}

impl<F: Field + TwoAdicField + PrimeField32> PartialEq for PolynomialCoeffs<F> {
    fn eq(&self, other: &Self) -> bool {
        let max_terms = self.coeffs.len().max(other.coeffs.len());
        for i in 0..max_terms {
            let self_i = self.coeffs.get(i).cloned().unwrap_or(F::ZERO);
            let other_i = other.coeffs.get(i).cloned().unwrap_or(F::ZERO);
            if self_i != other_i {
                return false;
            }
        }
        true
    }
}

impl<F: Field + TwoAdicField + PrimeField32> Eq for PolynomialCoeffs<F> {}

impl<F: Field + TwoAdicField + PrimeField32> From<Vec<F>> for PolynomialCoeffs<F> {
    fn from(coeffs: Vec<F>) -> Self {
        Self::new(coeffs)
    }
}

impl<F: Field + TwoAdicField + PrimeField32> Add for &PolynomialCoeffs<F> {
    type Output = PolynomialCoeffs<F>;

    fn add(self, rhs: Self) -> Self::Output {
        let len = max(self.len(), rhs.len());
        let a = self.padded(len).coeffs;
        let b = rhs.padded(len).coeffs;
        let coeffs = a.into_iter().zip(b).map(|(x, y)| x + y).collect();
        PolynomialCoeffs::new(coeffs)
    }
}

impl<F: Field + TwoAdicField + PrimeField32> Sum for PolynomialCoeffs<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::empty(), |acc, p| &acc + &p)
    }
}

impl<F: Field + TwoAdicField + PrimeField32> Sub for &PolynomialCoeffs<F> {
    type Output = PolynomialCoeffs<F>;

    fn sub(self, rhs: Self) -> Self::Output {
        let len = max(self.len(), rhs.len());
        let mut coeffs = self.padded(len).coeffs;
        for (i, &c) in rhs.coeffs.iter().enumerate() {
            coeffs[i] -= c;
        }
        PolynomialCoeffs::new(coeffs)
    }
}

impl<F: Field + TwoAdicField + PrimeField32> AddAssign for PolynomialCoeffs<F> {
    fn add_assign(&mut self, rhs: Self) {
        let len = max(self.len(), rhs.len());
        self.coeffs.resize(len, F::ZERO);
        for (l, r) in self.coeffs.iter_mut().zip(rhs.coeffs) {
            *l += r;
        }
    }
}

impl<F: Field + TwoAdicField + PrimeField32> AddAssign<&Self> for PolynomialCoeffs<F> {
    fn add_assign(&mut self, rhs: &Self) {
        let len = max(self.len(), rhs.len());
        self.coeffs.resize(len, F::ZERO);
        for (l, &r) in self.coeffs.iter_mut().zip(&rhs.coeffs) {
            *l += r;
        }
    }
}

impl<F: Field + TwoAdicField + PrimeField32> SubAssign for PolynomialCoeffs<F> {
    fn sub_assign(&mut self, rhs: Self) {
        let len = max(self.len(), rhs.len());
        self.coeffs.resize(len, F::ZERO);
        for (l, r) in self.coeffs.iter_mut().zip(rhs.coeffs) {
            *l -= r;
        }
    }
}

impl<F: Field + TwoAdicField + PrimeField32> SubAssign<&Self> for PolynomialCoeffs<F> {
    fn sub_assign(&mut self, rhs: &Self) {
        let len = max(self.len(), rhs.len());
        self.coeffs.resize(len, F::ZERO);
        for (l, &r) in self.coeffs.iter_mut().zip(&rhs.coeffs) {
            *l -= r;
        }
    }
}

impl<F: Field + TwoAdicField + PrimeField32> Mul<F> for &PolynomialCoeffs<F> {
    type Output = PolynomialCoeffs<F>;

    fn mul(self, rhs: F) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|&x| rhs * x).collect();
        PolynomialCoeffs::new(coeffs)
    }
}

impl<F: Field + TwoAdicField + PrimeField32> MulAssign<F> for PolynomialCoeffs<F> {
    fn mul_assign(&mut self, rhs: F) {
        self.coeffs.iter_mut().for_each(|x| *x *= rhs);
    }
}

impl<F: Field + TwoAdicField + PrimeField32> Mul for &PolynomialCoeffs<F> {
    type Output = PolynomialCoeffs<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: Self) -> Self::Output {
        let new_len = (self.len() + rhs.len()).next_power_of_two();
        let a = self.padded(new_len);
        let b = rhs.padded(new_len);
        let a_evals = a.fft();
        let b_evals = b.fft();

        let mul_evals: Vec<F> = a_evals
            .values
            .into_iter()
            .zip(b_evals.values)
            .map(|(pa, pb)| pa * pb)
            .collect();
        ifft(mul_evals.into())
    }
}

