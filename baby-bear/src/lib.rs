//! The prime field `F_p` where `p = 2^31 - 2^27 + 1`.

#![no_std]


use core::fmt;
use core::fmt::{Debug, Display, Formatter};
use core::hash::{Hash, Hasher};
use core::iter::{Product, Sum};
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use p3_field::{AbstractField, Field, PrimeField, PrimeField32};
use rand::distributions::{Distribution, Standard};
use rand::Rng;

/// The prime field `F_p` where `p = 2^31 - 2^27 + 1`.
#[derive(Copy, Clone, Default)]
pub struct BabyBear {
    /// Not necessarily canonical, but must fit in 31 bits.
    value: u32,
}

impl BabyBear {
    const fn new(value: u32) -> Self {
        Self { value: encode(value % P) }
    }

    pub fn new_raw(value: u32) -> Self {
        Self { value: value }
    }


    /// Raise to a power of `n`.
    fn pow(self, n: usize) -> Self {
        let mut n = n;
        let mut tot = Self::ONE;
        let mut x = self.clone();
        while n != 0 {
            if n % 2 == 1 {
                tot *= x;
            }
            n = n / 2;
            x *= x;
        }
        tot
    }
}

impl PartialEq for BabyBear {
    fn eq(&self, other: &Self) -> bool {
        self.as_canonical_u32() == other.as_canonical_u32()
    }
}

impl Eq for BabyBear {}

impl Hash for BabyBear {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u32(self.as_canonical_u32());
    }
}

impl Ord for BabyBear {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.as_canonical_u32().cmp(&other.as_canonical_u32())
    }
}

impl PartialOrd for BabyBear {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for BabyBear {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Display::fmt(&decode(self.value), f)
    }
}

impl Debug for BabyBear {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Debug::fmt(&decode(self.value), f)
    }
}

impl Distribution<BabyBear> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BabyBear {
        loop {
            let next_u31 = rng.next_u32() >> 1;
            let is_canonical = next_u31 < BabyBear::ORDER_U32 << 1;
            if is_canonical {
                return BabyBear::new(next_u31 % BabyBear::ORDER_U32);
            }
        }
    }
}

/// The modulus of the field.
const P: u32 = 15 * (1 << 27) + 1;

/// The modulus of the field as a u64.
const P_U64: u64 = P as u64;

/// Wrapping addition of [Elem] using Baby Bear field modulus
fn add(lhs: u32, rhs: u32) -> u32 {
    let x = lhs.wrapping_add(rhs);
    return if x >= P { x - P } else { x };
}

/// Wrapping subtraction of [Elem] using Baby Bear field modulus
fn sub(lhs: u32, rhs: u32) -> u32 {
    let x = lhs.wrapping_sub(rhs);
    return if x > P { x.wrapping_add(P) } else { x };
}

/// Wrapping multiplication of [Elem]  using Baby Bear field modulus
// Copied from the C++ implementation (fp.h)
const fn mul(lhs: u32, rhs: u32) -> u32 {
    // uint64_t o64 = uint64_t(a) * uint64_t(b);
    let mut o64: u64 = (lhs as u64).wrapping_mul(rhs as u64);
    // uint32_t low = -uint32_t(o64);
    let low: u32 = 0u32.wrapping_sub(o64 as u32);
    // uint32_t red = M * low;
    let red = M.wrapping_mul(low);
    // o64 += uint64_t(red) * uint64_t(P);
    o64 += (red as u64).wrapping_mul(P_U64);
    // uint32_t ret = o64 >> 32;
    let ret = (o64 >> 32) as u32;
    // return (ret >= P ? ret - P : ret);
    if ret >= P {
        ret - P
    } else {
        ret
    }
}

pub fn mul_for_test_1000(lhs: u32, rhs: u32) -> u32 {
    let mut ret = lhs;
    for _ in 0..1000 {
        ret = mul(ret, rhs);
    }
    ret
}

pub fn mul_with_two_muls_1000(a: u32, b: u32) -> u32 {
    let mut ret = a;
    for _ in 0..1000 {
        // uint64_t o64 = uint64_t(a) * uint64_t(b);
        let mut o64: u64 = (ret as u64).wrapping_mul(b as u64);
        // uint32_t low = -uint32_t(o64);
        let _low: u32 = 0u32.wrapping_sub(o64 as u32);
        // uint32_t red = M * low;
        // let red = M; //.wrapping_mul(low);
        // o64 += uint64_t(red) * uint64_t(P);
        o64 += (M as u64).wrapping_mul(P_U64);
        // uint32_t ret = o64 >> 32;
        ret = (o64 >> 32) as u32;
        // return (ret >= P ? ret - P : ret);
        if ret >= P {
            ret = ret - P;
        }
    }
    ret
}

pub fn u32_mul_1000(lhs: u32, rhs: u32) -> u32 {
    let mut ret = lhs;
    for _ in 0..1000 {
        ret = ret.wrapping_mul(rhs) + 1;
    }
    ret
}

pub fn u32_in_64_mul_1000(lhs: u32, rhs: u32) -> u64 {
    let mut ret = lhs;
    for _ in 0..1000 {
        ret = ((ret as u64).wrapping_mul(rhs as u64) + 1) as u32;
    }
    ret as u64
}

pub fn u64_mul_1000(lhs: u64, rhs: u64) -> u64 {
    let mut ret = lhs;
    for _ in 0..1000 {
        ret = ret.wrapping_mul(rhs) + 1;
    }
    ret
}

const M: u32 = 0x88000001;
const R2: u32 = 1172168163;

/// Encode to Montgomery form from direct form.
const fn encode(a: u32) -> u32 {
    mul(R2, a)
}

/// Decode from Montgomery form from direct form.
const fn decode(a: u32) -> u32 {
    mul(1, a)
}

impl AbstractField for BabyBear {
    const ZERO: Self = Self::new(0);
    const ONE: Self = Self::new(1);
    const TWO: Self = Self::new(2);
    const NEG_ONE: Self = Self::new(Self::ORDER_U32 - 1);

    fn from_canonical_u8(n: u8) -> Self {
        Self::new(u32::from(n))
    }

    fn from_canonical_u16(n: u8) -> Self {
        Self::new(u32::from(n))
    }

    fn from_canonical_u32(n: u32) -> Self {
        Self::new(n)
    }

    /// Convert from `u64`. Undefined behavior if the input is outside the canonical range.
    fn from_canonical_u64(n: u64) -> Self {
        Self::new(
            n.try_into()
                .expect("Too large to be a canonical BabyBear encoding"),
        )
    }

    /// Convert from `usize`. Undefined behavior if the input is outside the canonical range.
    fn from_canonical_usize(n: usize) -> Self {
        Self::new(
            n.try_into()
                .expect("Too large to be a canonical BabyBear encoding"),
        )
    }

    fn from_wrapped_u32(n: u32) -> Self {
        Self {
            value: encode(n)
        }
    }

    fn from_wrapped_u64(_n: u64) -> Self {
        todo!()
    }

    // Sage: GF(2^31 - 2^27 + 1).multiplicative_generator()
    fn multiplicative_group_generator() -> Self {
        Self::new(137)
    }
}

impl Field for BabyBear {
    // TODO: Add cfg-guarded Packing for AVX2, NEON, etc.
    type Packing = Self;

    fn is_zero(&self) -> bool {
        self.value == 0 || self.value == Self::ORDER_U32 || self.value == Self::ORDER_U32 * 2
    }

    fn mul_2exp_u64(&self, exp: u64) -> Self {
        let exp = (exp % 31) as u8;
        *self * Self::from_canonical_u32(1 << exp)
    }

    fn div_2exp_u64(&self, exp: u64) -> Self {
        let exp = (exp % 31) as u8;
        *self * Self::from_canonical_u32(1 << exp).try_inverse().unwrap()
    }

    fn try_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        Some(self.pow((P - 2) as usize))
    }
}

impl PrimeField for BabyBear {}

impl PrimeField32 for BabyBear {
    const ORDER_U32: u32 = (1 << 31) - (1 << 27) + 1;
    const CHARACTERISTIC_TWO_ADICITY: u32 = 27;

    fn as_canonical_u32(&self) -> u32 {
        if self.value >= Self::ORDER_U32 {
            self.value % Self::ORDER_U32
        } else {
            self.value
        }
    }
}

impl Add for BabyBear {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            value: add(self.value, rhs.value)
        }
    }
}

impl AddAssign for BabyBear {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sum for BabyBear {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|x, y| x + y).unwrap_or(Self::ZERO)
    }
}

impl Sub for BabyBear {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            value: sub(self.value, rhs.value)
        }
    }
}

impl SubAssign for BabyBear {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for BabyBear {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(Self::ORDER_U32 - (self.value % Self::ORDER_U32))
    }
}

impl Mul for BabyBear {
    type Output = Self;

    #[allow(clippy::cast_possible_truncation)]
    fn mul(self, rhs: Self) -> Self {
        Self { value: mul(self.value, rhs.value) }
    }
}

impl MulAssign for BabyBear {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Product for BabyBear {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|x, y| x * y).unwrap_or(Self::ONE)
    }
}

impl Div for BabyBear {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self {
        self * rhs.try_inverse().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use crate::BabyBear;
    use p3_field::{AbstractField, Field, PrimeField32};

    type F = BabyBear;

    // These tests are for Mersenne fields, not baby bear
    // So they will not pass for now.

    #[test]
    fn add() {
        assert_eq!(F::ONE + F::ONE, F::TWO);
        assert_eq!(F::NEG_ONE + F::ONE, F::ZERO);
        assert_eq!(F::NEG_ONE + F::TWO, F::ONE);
        assert_eq!(F::NEG_ONE + F::NEG_ONE, F::new(F::ORDER_U32 - 2));
    }

    #[test]
    fn sub() {
        assert_eq!(F::ONE - F::ONE, F::ZERO);
        assert_eq!(F::TWO - F::TWO, F::ZERO);
        assert_eq!(F::NEG_ONE - F::NEG_ONE, F::ZERO);
        assert_eq!(F::TWO - F::ONE, F::ONE);
        assert_eq!(F::NEG_ONE - F::ZERO, F::NEG_ONE);
    }

    #[test]
    fn mul_2exp_u64() {
        // 1 * 2^0 = 1.
        assert_eq!(F::ONE.mul_2exp_u64(0), F::ONE);
        // 2 * 2^30 = 2^31 = 1.
        assert_eq!(F::TWO.mul_2exp_u64(30), F::ONE);
        // 5 * 2^2 = 20.
        assert_eq!(F::new(5).mul_2exp_u64(2), F::new(20));
    }

    #[test]
    fn div_2exp_u64() {
        // 1 / 2^0 = 1.
        assert_eq!(F::ONE.div_2exp_u64(0), F::ONE);
        // 2 / 2^0 = 2.
        assert_eq!(F::TWO.div_2exp_u64(0), F::TWO);
        // 32 / 2^5 = 1.
        assert_eq!(F::new(32).div_2exp_u64(5), F::new(1));
    }

    #[test]
    fn inverse() {
        assert_eq!(F::new(172).inverse() * F::new(172), F::ONE);
    }
}
