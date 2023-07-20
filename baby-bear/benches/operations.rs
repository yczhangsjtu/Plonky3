use p3_field::PrimeField32;
use p3_field::AbstractField;
use core::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use p3_baby_bear::BabyBear;
use p3_field::Field;
use rand::distributions::{Distribution, Standard};

use std::any::type_name;


fn bench_baby_bear_operations(c: &mut Criterion) {
    baby_bear_operations::<BabyBear>(c);
}

fn baby_bear_operations<F: Field>(c: &mut Criterion)
where
    Standard: Distribution<F>,
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<F>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    const BATCH_SIZE: u32 = 1 << 12;

    let a = F::from_canonical_u32(1234567890);
    let b = F::from_canonical_u32(2345678901);

    group.bench_with_input(BenchmarkId::from_parameter("baby mul 1234567890"), &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let mut a = *a;
            let b = *b;
            for _ in 0..BATCH_SIZE {
                a = a * b;
            }
            a
        });
    });

    let a = F::from_canonical_u32(1567890123);
    let b = F::from_canonical_u32(3890123456);

    group.bench_with_input(BenchmarkId::from_parameter("baby mul 1567890123"), &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let mut a = *a;
            let b = *b;
            for _ in 0..BATCH_SIZE {
                a = a * b;
            }
            a
        });
    });

    let a: u32 = BabyBear::from_canonical_u32(1234567890).as_canonical_u32();
    let b: u32 = BabyBear::from_canonical_u32(2345678901).as_canonical_u32();

    group.bench_with_input(BenchmarkId::from_parameter("u32 in u64 mul"), &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let mut a = *a;
            let b = *b;
            for _ in 0..BATCH_SIZE {
                a = (((a as u64).wrapping_mul(b as u64) >> 32) as u32) + 1234567;
            }
            a
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("u64 mul"), &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let mut a = (*a as u64) * (*b as u64);
            let b = a;
            for _ in 0..BATCH_SIZE {
                a = a * b + 123456789012345;
            }
            a
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby direct"), &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let mut a = *a;
            let b = *b;
            const M: u32 = 0x88000001;
            const P: u32 = 15 * (1 << 27) + 1;
            const P_U64: u64 = P as u64;
            for _ in 0..BATCH_SIZE {
                // uint64_t o64 = uint64_t(a) * uint64_t(b);
                let mut o64: u64 = (a as u64).wrapping_mul(b as u64);
                // uint32_t low = -uint32_t(o64);
                let low: u32 = 0u32.wrapping_sub(o64 as u32);
                // uint32_t red = M * low;
                let red = M.wrapping_mul(low);
                // o64 += uint64_t(red) * uint64_t(P);
                o64 += (red as u64).wrapping_mul(P_U64);
                // uint32_t ret = o64 >> 32;
                a = (o64 >> 32) as u32;
                if a >= P {
                    a = a - P;
                }
            }
            a
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby direct 2 muls"), &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let mut a = *a;
            let b = *b;
            const M: u32 = 0x88000001;
            const P: u32 = 15 * (1 << 27) + 1;
            const P_U64: u64 = P as u64;
            for _ in 0..BATCH_SIZE {
                // uint64_t o64 = uint64_t(a) * uint64_t(b);
                let mut o64: u64 = ((a as u64) << 32) + (b as u64);
                // uint32_t low = -uint32_t(o64);
                let low: u32 = 0u32.wrapping_sub(o64 as u32);
                // uint32_t red = M * low;
                let red = M.wrapping_mul(low);
                // o64 += uint64_t(red) * uint64_t(P);
                o64 += (red as u64).wrapping_mul(P_U64);
                // uint32_t ret = o64 >> 32;
                a = (o64 >> 32) as u32;
                if a >= P {
                    a = a - P;
                }
            }
            a
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby inverse"), &a, |bench, a| {
        bench.iter(|| {
            let mut a = BabyBear::new_raw(*a);
            for _ in 0..BATCH_SIZE {
                a = a.try_inverse().unwrap();
            }
            a
        });
    });
}

criterion_group!(benches, bench_baby_bear_operations);
criterion_main!(benches);
