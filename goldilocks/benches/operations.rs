use rand::Rng;
use p3_field::PrimeField64;
use core::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use p3_goldilocks::Goldilocks;
use p3_field::Field;
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;

use std::any::type_name;


fn bench_goldilocks_operations(c: &mut Criterion) {
    goldilocks_operations::<Goldilocks>(c);
}

fn goldilocks_operations<F: Field + PrimeField64>(c: &mut Criterion)
where
    Standard: Distribution<F>,
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<F>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    let mut rng = thread_rng();

    let a = F::from_canonical_u64(12345678902345678901);
    let b = F::from_canonical_u64(13456789011234567890);

    let mut av = vec![F::ZERO; 1000];
    let mut bv = vec![F::ZERO; 1000];
    for i in 0..1000 {
        av[i] = F::from_wrapped_u64(rng.gen::<u64>());
        bv[i] = F::from_wrapped_u64(rng.gen::<u64>());
    }

    group.bench_with_input(BenchmarkId::from_parameter("goldilocks mul vec 1234567890"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1000];
        bench.iter(|| {
            for i in 0..1000 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("goldilocks add vec"), &(av, bv), |bench, (av, bv)| {
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
        av[i] = F::from_wrapped_u64(rng.gen::<u64>());
        bv[i] = F::from_wrapped_u64(rng.gen::<u64>());
    }

    group.bench_with_input(BenchmarkId::from_parameter("goldilocks mul vec 1 << 14"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
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
        av[i] = F::from_wrapped_u64(rng.gen::<u64>());
        bv[i] = F::from_wrapped_u64(rng.gen::<u64>());
    }

    group.bench_with_input(BenchmarkId::from_parameter("goldilocks mul vec 1 << 16"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1 << 16];
        bench.iter(|| {
            for i in 0..1 << 16 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("goldilocks mul 12345678902345678901"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a * b;
        });
    });

    let a = F::from_canonical_u64(15678901233890123456);
    let b = F::from_canonical_u64(12901234561567890123);

    group.bench_with_input(BenchmarkId::from_parameter("goldilocks mul 15678901233890123456"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a * b;
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("goldilocks add"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a + b;
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("u64 in u128 mul"), &(a, b), |bench, (a, b)| {
        let mut a = a.as_canonical_u64();
        let b = b.as_canonical_u64();
        bench.iter(|| {
            a = ((((a as u128) * (b as u128)) >> 64) as u64) + 123456789012345;
        });
    });
    
    group.bench_with_input(BenchmarkId::from_parameter("goldilocks inverse"), &a, |bench, a| {
        let mut a = *a;
        bench.iter(|| {
            a = a.try_inverse().unwrap();
        });
    });
}

criterion_group!(benches, bench_goldilocks_operations);
criterion_main!(benches);
