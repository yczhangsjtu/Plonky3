use rand::Rng;
use core::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use p3_mersenne_31::Mersenne31;
use p3_field::Field;
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;

use std::any::type_name;


fn bench_baby_bear_operations(c: &mut Criterion) {
    baby_bear_operations::<Mersenne31>(c);
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

    group.bench_with_input(BenchmarkId::from_parameter("mersenne mul vec 1234567890"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1000];
        bench.iter(|| {
            for i in 0..1000 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("mersenne add vec"), &(av, bv), |bench, (av, bv)| {
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

    group.bench_with_input(BenchmarkId::from_parameter("mersenne mul vec 1 << 14"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
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

    group.bench_with_input(BenchmarkId::from_parameter("mersenne mul vec 1 << 16"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1 << 16];
        bench.iter(|| {
            for i in 0..1 << 16 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("mersenne 1234567890"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a * b;
        });
    });

    let a = F::from_canonical_u32(1567890123);
    let b = F::from_canonical_u32(3890123456);

    group.bench_with_input(BenchmarkId::from_parameter("mersenne 1567890123"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a * b;
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("mersenne add"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a + b;
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("mersenne inverse"), &a, |bench, a| {
        let mut a = *a;
        bench.iter(|| {
            a = a.try_inverse().unwrap();
        });
    });
}

criterion_group!(benches, bench_baby_bear_operations);
criterion_main!(benches);
