use core::time::Duration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use p3_brakedown::fast_registry;
use p3_code::CodeOrFamily;
use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;
use p3_goldilocks::Goldilocks;
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;
use std::any::type_name;

const BATCH_SIZE: usize = 1 << 12;

fn bench_encode_goldilocks(c: &mut Criterion) {
    encode::<Goldilocks, 20>(c);
}

fn encode<F: Field, const ROW_WEIGHT: usize>(c: &mut Criterion)
where
    Standard: Distribution<F>,
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<F>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    let mut rng = thread_rng();
    for n_log in [14, 16] {
        let n = 1 << n_log;

        let code = fast_registry();

        let mut messages = RowMajorMatrix::rand(&mut rng, n, BATCH_SIZE);

        group.bench_with_input(BenchmarkId::from_parameter(n), &code, |b, code| {
            b.iter(|| {
                messages.values.truncate(n * BATCH_SIZE);
                code.encode_batch(messages.clone());
            });
        });
    }
}

criterion_group!(benches, bench_encode_goldilocks);
criterion_main!(benches);
