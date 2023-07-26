use p3_mersenne_31::Mersenne31Complex;
use rand::Rng;
use core::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use p3_mersenne_31::Mersenne31;
use p3_field::Field;
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;

use std::any::type_name;


fn bench_mersenne_complex_operations(c: &mut Criterion) {
    mersenne_complex_operations::<Mersenne31Complex<Mersenne31>>(c);
}

fn mersenne_complex_operations<F: Field>(c: &mut Criterion)
where
    Standard: Distribution<F>,
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<F>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    let mut rng = thread_rng();

    let a: F = rng.gen();
    let b: F = rng.gen();

    let mut av = vec![F::ZERO; 1000];
    let mut bv = vec![F::ZERO; 1000];
    for i in 0..1000 {
        av[i] = rng.gen();
        bv[i] = rng.gen();
    }

    group.bench_with_input(BenchmarkId::from_parameter("complex mul vec 1234567890"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1000];
        bench.iter(|| {
            for i in 0..1000 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("complex add vec"), &(av, bv), |bench, (av, bv)| {
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
        av[i] = rng.gen();
        bv[i] = rng.gen();
    }

    group.bench_with_input(BenchmarkId::from_parameter("complex mul vec 1 << 14"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
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
        av[i] = rng.gen();
        bv[i] = rng.gen();
    }

    group.bench_with_input(BenchmarkId::from_parameter("complex mul vec 1 << 16"), &(av.clone(), bv.clone()), |bench, (av, bv)| {
        let mut v = vec![F::ZERO; 1 << 16];
        bench.iter(|| {
            for i in 0..1 << 16 {
                v[i] = av[i] * bv[i];
            }
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("complex mul"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a * b;
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("complex add"), &(a, b), |bench, (a, b)| {
        let mut a = *a;
        let b = *b;
        bench.iter(|| {
            a = a + b;
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("complex inverse"), &a, |bench, a| {
        let mut a = *a;
        bench.iter(|| {
            a = a.try_inverse().unwrap();
        });
    });
}

criterion_group!(benches, bench_mersenne_complex_operations);
criterion_main!(benches);
