//! The prime field known as Goldilocks, defined as `F_p` where `p = 2^64 - 2^32 + 1`.

#![no_std]

use core::fmt;
use core::fmt::{Debug, Display, Formatter};
use core::hash::{Hash, Hasher};
use core::iter::{Product, Sum};
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use p3_field::{AbstractField, Field, PrimeField, PrimeField64, TwoAdicField};
use p3_util::{assume, branch_hint};
use rand::distributions::{Distribution, Standard};
use rand::Rng;

mod inversion;
use inversion::try_inverse_u64;

/// The prime field known as Goldilocks, defined as `F_p` where `p = 2^64 - 2^32 + 1`.
#[derive(Copy, Clone, Default)]
pub struct Goldilocks {
    /// Not necessarily canonical.
    value: u64,
}

impl Goldilocks {
    const fn new(value: u64) -> Self {
        Self { value }
    }

    /// Two's complement of `ORDER`, i.e. `2^64 - ORDER = 2^32 - 1`.
    const NEG_ORDER: u64 = Self::ORDER_U64.wrapping_neg();
}

impl PartialEq for Goldilocks {
    fn eq(&self, other: &Self) -> bool {
        self.as_canonical_u64() == other.as_canonical_u64()
    }
}

impl Eq for Goldilocks {}

impl Hash for Goldilocks {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.as_canonical_u64());
    }
}

impl Ord for Goldilocks {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.as_canonical_u64().cmp(&other.as_canonical_u64())
    }
}

impl PartialOrd for Goldilocks {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for Goldilocks {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Display::fmt(&self.value, f)
    }
}

impl Debug for Goldilocks {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Debug::fmt(&self.value, f)
    }
}

impl Distribution<Goldilocks> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Goldilocks {
        loop {
            let next_u64 = rng.next_u64();
            let is_canonical = next_u64 < Goldilocks::ORDER_U64;
            if is_canonical {
                return Goldilocks::new(next_u64);
            }
        }
    }
}

impl AbstractField for Goldilocks {
    const ZERO: Self = Self::new(0);
    const ONE: Self = Self::new(1);
    const TWO: Self = Self::new(2);
    const NEG_ONE: Self = Self::new(Self::ORDER_U64 - 1);

    fn from_bool(b: bool) -> Self {
        Self::new(u64::from(b))
    }

    fn from_canonical_u8(n: u8) -> Self {
        Self::new(u64::from(n))
    }

    fn from_canonical_u16(n: u16) -> Self {
        Self::new(u64::from(n))
    }

    fn from_canonical_u32(n: u32) -> Self {
        Self::new(u64::from(n))
    }

    fn from_canonical_u64(n: u64) -> Self {
        Self::new(n)
    }

    fn from_canonical_usize(n: usize) -> Self {
        Self::new(n as u64)
    }

    fn from_wrapped_u32(n: u32) -> Self {
        // A u32 must be canonical, plus we don't store canonical encodings anyway, so there's no
        // need for a reduction.
        Self::new(u64::from(n))
    }

    fn from_wrapped_u64(n: u64) -> Self {
        // There's no need to reduce `n` to canonical form, as our internal encoding is
        // non-canonical, so there's no need for a reduction.
        Self::new(n)
    }

    // Sage: GF(2^64 - 2^32 + 1).multiplicative_generator()
    fn multiplicative_group_generator() -> Self {
        Self::new(7)
    }
}

impl Field for Goldilocks {
    // TODO: Add cfg-guarded Packing for AVX2, NEON, etc.
    type Packing = Self;

    fn is_zero(&self) -> bool {
        self.value == 0 || self.value == Self::ORDER_U64
    }

    fn try_inverse(&self) -> Option<Self> {
        try_inverse_u64(self)
    }
}

impl PrimeField for Goldilocks {}

impl PrimeField64 for Goldilocks {
    const ORDER_U64: u64 = 0xFFFF_FFFF_0000_0001;
    const CHARACTERISTIC_TWO_ADICITY: u64 = 32;

    fn as_canonical_u64(&self) -> u64 {
        let mut c = self.value;
        // We only need one condition subtraction, since 2 * ORDER would not fit in a u64.
        if c >= Self::ORDER_U64 {
            c -= Self::ORDER_U64;
        }
        c
    }
}

impl TwoAdicField for Goldilocks {
    const TWO_ADICITY: usize = 32;

    fn power_of_two_generator() -> Self {
        Self::new(1_753_635_133_440_165_772)
    }

    /// Compute the inverse of 2^exp in this field.
    #[inline]
    fn inverse_2exp(exp: usize) -> Self {
        // Let p = char(F). Since 2^exp is in the prime subfield, i.e. an
        // element of GF_p, its inverse must be as well. Thus we may add
        // multiples of p without changing the result. In particular,
        // 2^-exp = 2^-exp - p 2^-exp
        //        = 2^-exp (1 - p)
        //        = p - (p - 1) / 2^exp

        // If this field's two adicity, t, is at least exp, then 2^exp divides
        // p - 1, so this division can be done with a simple bit shift. If
        // exp > t, we repeatedly multiply by 2^-t and reduce exp until it's in
        // the right range.

        let p = Self::ORDER_U64;
        // NB: The only reason this is split into two cases is to save
        // the multiplication (and possible calculation of
        // inverse_2_pow_adicity) in the usual case that exp <=
        // TWO_ADICITY. Can remove the branch and simplify if that
        // saving isn't worth it.

        if exp > Self::CHARACTERISTIC_TWO_ADICITY as usize {
            // NB: This should be a compile-time constant
            let inverse_2_pow_adicity: Self =
                Self::from_canonical_u64(p - ((p - 1) >> Self::CHARACTERISTIC_TWO_ADICITY));

            let mut res = inverse_2_pow_adicity;
            let mut e = exp - Self::CHARACTERISTIC_TWO_ADICITY as usize;

            while e > Self::CHARACTERISTIC_TWO_ADICITY as usize {
                res *= inverse_2_pow_adicity;
                e -= Self::CHARACTERISTIC_TWO_ADICITY as usize;
            }
            res * Self::from_canonical_u64(p - ((p - 1) >> e))
        } else {
            Self::from_canonical_u64(p - ((p - 1) >> exp))
        }
    }
}

impl Add for Goldilocks {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let (sum, over) = self.value.overflowing_add(rhs.value);
        let (mut sum, over) = sum.overflowing_add(u64::from(over) * Self::NEG_ORDER);
        if over {
            // NB: self.value > Self::ORDER && rhs.value > Self::ORDER is necessary but not
            // sufficient for double-overflow.
            // This assume does two things:
            //  1. If compiler knows that either self.value or rhs.value <= ORDER, then it can skip
            //     this check.
            //  2. Hints to the compiler how rare this double-overflow is (thus handled better with
            //     a branch).
            assume(self.value > Self::ORDER_U64 && rhs.value > Self::ORDER_U64);
            branch_hint();
            sum += Self::NEG_ORDER; // Cannot overflow.
        }
        Self::new(sum)
    }
}

impl AddAssign for Goldilocks {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sum for Goldilocks {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|x, y| x + y).unwrap_or(Self::ZERO)
    }
}

impl Sub for Goldilocks {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let (diff, under) = self.value.overflowing_sub(rhs.value);
        let (mut diff, under) = diff.overflowing_sub(u64::from(under) * Self::NEG_ORDER);
        if under {
            // NB: self.value < NEG_ORDER - 1 && rhs.value > ORDER is necessary but not
            // sufficient for double-underflow.
            // This assume does two things:
            //  1. If compiler knows that either self.value >= NEG_ORDER - 1 or rhs.value <= ORDER,
            //     then it can skip this check.
            //  2. Hints to the compiler how rare this double-underflow is (thus handled better
            //     with a branch).
            assume(self.value < Self::NEG_ORDER - 1 && rhs.value > Self::ORDER_U64);
            branch_hint();
            diff -= Self::NEG_ORDER; // Cannot underflow.
        }
        Self::new(diff)
    }
}

impl SubAssign for Goldilocks {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for Goldilocks {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(Self::ORDER_U64 - self.as_canonical_u64())
    }
}

impl Mul for Goldilocks {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        reduce128(u128::from(self.value) * u128::from(rhs.value))
    }
}

impl MulAssign for Goldilocks {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Product for Goldilocks {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|x, y| x * y).unwrap_or(Self::ONE)
    }
}

impl Div for Goldilocks {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

/// Reduces to a 64-bit value. The result might not be in canonical form; it could be in between the
/// field order and `2^64`.
#[inline]
fn reduce128(x: u128) -> Goldilocks {
    let (x_lo, x_hi) = split(x); // This is a no-op
    let x_hi_hi = x_hi >> 32;
    let x_hi_lo = x_hi & Goldilocks::NEG_ORDER;

    let (mut t0, borrow) = x_lo.overflowing_sub(x_hi_hi);
    if borrow {
        branch_hint(); // A borrow is exceedingly rare. It is faster to branch.
        t0 -= Goldilocks::NEG_ORDER; // Cannot underflow.
    }
    let t1 = x_hi_lo * Goldilocks::NEG_ORDER;
    let t2 = unsafe { add_no_canonicalize_trashing_input(t0, t1) };
    Goldilocks::new(t2)
}

#[inline]
#[allow(clippy::cast_possible_truncation)]
fn split(x: u128) -> (u64, u64) {
    (x as u64, (x >> 64) as u64)
}

/// Fast addition modulo ORDER for x86-64.
/// This function is marked unsafe for the following reasons:
///   - It is only correct if x + y < 2**64 + ORDER = 0x1ffffffff00000001.
///   - It is only faster in some circumstances. In particular, on x86 it overwrites both inputs in
///     the registers, so its use is not recommended when either input will be used again.
#[inline(always)]
#[cfg(target_arch = "x86_64")]
unsafe fn add_no_canonicalize_trashing_input(x: u64, y: u64) -> u64 {
    let res_wrapped: u64;
    let adjustment: u64;
    core::arch::asm!(
        "add {0}, {1}",
        // Trick. The carry flag is set iff the addition overflowed.
        // sbb x, y does x := x - y - CF. In our case, x and y are both {1:e}, so it simply does
        // {1:e} := 0xffffffff on overflow and {1:e} := 0 otherwise. {1:e} is the low 32 bits of
        // {1}; the high 32-bits are zeroed on write. In the end, we end up with 0xffffffff in {1}
        // on overflow; this happens be NEG_ORDER.
        // Note that the CPU does not realize that the result of sbb x, x does not actually depend
        // on x. We must write the result to a register that we know to be ready. We have a
        // dependency on {1} anyway, so let's use it.
        "sbb {1:e}, {1:e}",
        inlateout(reg) x => res_wrapped,
        inlateout(reg) y => adjustment,
        options(pure, nomem, nostack),
    );
    assume(x != 0 || (res_wrapped == y && adjustment == 0));
    assume(y != 0 || (res_wrapped == x && adjustment == 0));
    // Add NEG_ORDER == subtract ORDER.
    // Cannot overflow unless the assumption if x + y < 2**64 + ORDER is incorrect.
    res_wrapped + adjustment
}

#[inline(always)]
#[cfg(not(target_arch = "x86_64"))]
unsafe fn add_no_canonicalize_trashing_input(x: u64, y: u64) -> u64 {
    let (res_wrapped, carry) = x.overflowing_add(y);
    // Below cannot overflow unless the assumption if x + y < 2**64 + ORDER is incorrect.
    res_wrapped + Goldilocks::NEG_ORDER * u64::from(carry)
}

#[cfg(test)]
mod tests {
    use super::*;

    type F = Goldilocks;

    #[test]
    fn test_goldilocks() {
        let f = F::new(100);
        assert_eq!(f.as_canonical_u64(), 100);

        // Over the Goldilocks field, the following set of equations hold
        // p               = 0
        // 2^64 - 2^32 + 1 = 0
        // 2^64            = 2^32 - 1
        let f = F::new(u64::MAX);
        assert_eq!(f.as_canonical_u64(), u32::MAX as u64 - 1);

        let f = F::from_canonical_u64(u64::MAX);
        assert_eq!(f.as_canonical_u64(), u32::MAX as u64 - 1);

        let f = F::from_canonical_u64(0);
        assert!(f.is_zero());

        let f = F::from_canonical_u64(F::ORDER_U64);
        assert!(f.is_zero());

        assert_eq!(
            F::multiplicative_group_generator().as_canonical_u64(),
            7_u64
        );

        let f_1 = F::new(1);
        let f_1_copy = F::new(1);

        let expected_result = F::ZERO;
        assert_eq!(f_1 - f_1_copy, expected_result);

        let expected_result = F::new(2);
        assert_eq!(f_1 + f_1_copy, expected_result);

        let f_2 = F::new(2);
        let expected_result = F::new(3);
        assert_eq!(f_1 + f_1_copy * f_2, expected_result);

        let expected_result = F::new(5);
        assert_eq!(f_1 + f_2 * f_2, expected_result);

        let f_p_minus_1 = F::from_canonical_u64(F::ORDER_U64 - 1);
        let expected_result = F::ZERO;
        assert_eq!(f_1 + f_p_minus_1, expected_result);

        let f_p_minus_2 = F::from_canonical_u64(F::ORDER_U64 - 2);
        let expected_result = F::from_canonical_u64(F::ORDER_U64 - 3);
        assert_eq!(f_p_minus_1 + f_p_minus_2, expected_result);

        let expected_result = F::new(1);
        assert_eq!(f_p_minus_1 - f_p_minus_2, expected_result);

        let expected_result = f_p_minus_1;
        assert_eq!(f_p_minus_2 - f_p_minus_1, expected_result);

        let expected_result = f_p_minus_2;
        assert_eq!(f_p_minus_1 - f_1, expected_result);

        let expected_result = F::new(3);
        assert_eq!(f_2 * f_2 - f_1, expected_result);

        // Generator check
        let expected_multiplicative_group_generator = F::new(7);
        assert_eq!(
            F::multiplicative_group_generator(),
            expected_multiplicative_group_generator
        );

        // Check on `reduce_u128`
        let x = u128::MAX;
        let y = reduce128(x);
        // The following equalitiy sequence holds, modulo p = 2^64 - 2^32 + 1
        // 2^128 - 1 = (2^64 - 1) * (2^64 + 1)
        //           = (2^32 - 1 - 1) * (2^32 - 1 + 1)
        //           = (2^32 - 2) * (2^32)
        //           = 2^64 - 2 * 2^32
        //           = 2^64 - 2^33
        //           = 2^32 - 1 - 2^33
        //           = - 2^32 - 1
        let expected_result = -F::new(2_u64.pow(32)) - F::new(1);
        assert_eq!(y, expected_result);
    }
}
