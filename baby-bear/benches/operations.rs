use rand::Rng;
use p3_field::PrimeField32;
use p3_field::AbstractField;
use core::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use p3_baby_bear::BabyBear;
use p3_baby_bear::{mul_for_test_1000, mul_with_two_muls_1000, u32_mul_1000, u64_mul_1000, u32_in_64_mul_1000};
use p3_field::Field;
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;

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
    let mut rng = thread_rng();

    let a = F::from_canonical_u32(1234567890);
    let b = F::from_canonical_u32(2345678901);
    let mut av = vec![F::ZERO; 1000];
    let mut bv = vec![F::ZERO; 1000];
    for i in 0..1000 {
        av[i] = F::from_wrapped_u32(rng.gen::<u32>());
        bv[i] = F::from_wrapped_u32(rng.gen::<u32>());
    }

    group.bench_with_input(BenchmarkId::from_parameter("baby mul vec 1234567890"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1000];
        bench.iter(|| {
            for i in 0..1000 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby add vec"), &(av, bv), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1000];
        bench.iter(|| {
            for i in 0..1000 {
                v[i] = av[i] + bv[i];
            }
        });
    });

    let mut av = vec![F::ZERO; 1 << 14];
    let mut bv = vec![F::ZERO; 1 << 14];
    for i in 0..1 << 14 {
        av[i] = F::from_wrapped_u32(rng.gen::<u32>());
        bv[i] = F::from_wrapped_u32(rng.gen::<u32>());
    }

    group.bench_with_input(BenchmarkId::from_parameter("baby mul vec 1 << 14"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1 << 14];
        bench.iter(|| {
            for i in 0..1 << 14 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    let mut av = vec![F::ZERO; 1 << 16];
    let mut bv = vec![F::ZERO; 1 << 16];
    for i in 0..1 << 16 {
        av[i] = F::from_wrapped_u32(rng.gen::<u32>());
        bv[i] = F::from_wrapped_u32(rng.gen::<u32>());
    }

    group.bench_with_input(BenchmarkId::from_parameter("baby mul vec 1 << 16"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1 << 16];
        bench.iter(|| {
            for i in 0..1 << 16 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby mul 1234567890"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a * b
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby add"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a + b
        });
    });

    let a: u32 = BabyBear::from_wrapped_u32(rng.gen::<u32>()).as_canonical_u32();
    let b: u32 = BabyBear::from_wrapped_u32(rng.gen::<u32>()).as_canonical_u32();

    group.bench_with_input(BenchmarkId::from_parameter("u32 in u64 mul"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = (u32_in_64_mul_1000(a, b) >> 32) as u32
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("u32 wrapping mul"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = u32_mul_1000(a, b)
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("u64 mul"), &(a, b), |bench, (a, b)| {
        let mut a = (*a as u64) * (*b as u64);
        let b = a;
        bench.iter(|| {
            a = u64_mul_1000(a, b)
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby direct"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = mul_for_test_1000(a, b)
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby direct 2 muls"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = mul_with_two_muls_1000(a, b)
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("baby inverse"), &a, |bench, a| {
        let a = BabyBear::new_raw(*a);
        bench.iter(|| {
            a.try_inverse().unwrap()
        });
    });
}

criterion_group!(benches, bench_baby_bear_operations);
criterion_main!(benches);
