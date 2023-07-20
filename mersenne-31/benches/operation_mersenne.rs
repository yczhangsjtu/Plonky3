use core::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use p3_mersenne_31::Mersenne31;
use p3_field::Field;
use rand::distributions::{Distribution, Standard};

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
    const BATCH_SIZE: u32 = 1 << 12;

    let a = F::from_canonical_u32(1234567890);
    let b = F::from_canonical_u32(2345678901);

    group.bench_with_input(BenchmarkId::from_parameter("mersenne 1234567890"), &(a, b), |bench, (a, b)| {
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

    group.bench_with_input(BenchmarkId::from_parameter("mersenne 1567890123"), &(a, b), |bench, (a, b)| {
        bench.iter(|| {
            let mut a = *a;
            let b = *b;
            for _ in 0..BATCH_SIZE {
                a = a * b;
            }
            a
        });
    });

    group.bench_with_input(BenchmarkId::from_parameter("mersenne inverse"), &a, |bench, a| {
        bench.iter(|| {
            let mut a = *a;
            for _ in 0..BATCH_SIZE {
                a = a.try_inverse().unwrap();
            }
            a
        });
    });
}

criterion_group!(benches, bench_baby_bear_operations);
criterion_main!(benches);
