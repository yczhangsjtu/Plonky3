use p3_field::AbstractField;
use p3_mersenne_31::Mersenne31Complex;
use p3_mersenne_31::Mersenne31;
use p3_field::PrimeField32;
use p3_field::PrimeField64;
use p3_goldilocks::Goldilocks;
use criterion::BenchmarkId;
use rand::Rng;
use p3_field::TwoAdicField;
use core::time::Duration;
use core::any::type_name;
use p3_field::Field;
use criterion::{criterion_group, criterion_main, Criterion};
use p3_baby_bear::BabyBear;
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;
use p3_lde::polynomial::PolynomialCoeffs;

fn bench_fft(c: &mut Criterion) {
    run_real_fft_complex(c);
    run_fft_complex::<Mersenne31Complex<Mersenne31>>(c);
    run_fft_64::<Goldilocks>(c);
    run_fft_32::<BabyBear>(c);
}

fn run_fft_32<F: Field + TwoAdicField + PrimeField32>(c: &mut Criterion)
where
    Standard: Distribution<F>,
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<F>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    let mut rng = thread_rng();

    let mut v = vec![F::ZERO; 1 << 14];
    for i in 0..(1 << 14) {
    	v[i] = F::from_wrapped_u32(rng.gen::<u32>());
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 14"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });

    let mut v = vec![F::ZERO; 1 << 16];
    for i in 0..(1 << 16) {
    	v[i] = F::from_wrapped_u32(rng.gen::<u32>());
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 16"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });

    let mut v = vec![F::ZERO; 1 << 20];
    for i in 0..(1 << 20) {
    	v[i] = F::from_wrapped_u32(rng.gen::<u32>());
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 20"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });
}


fn run_fft_64<F: Field + TwoAdicField + PrimeField64>(c: &mut Criterion)
where
    Standard: Distribution<F>,
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<F>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    let mut rng = thread_rng();

    let mut v = vec![F::ZERO; 1 << 14];
    for i in 0..(1 << 14) {
    	v[i] = F::from_wrapped_u64(rng.gen::<u64>());
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 14"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });

    let mut v = vec![F::ZERO; 1 << 16];
    for i in 0..(1 << 16) {
    	v[i] = F::from_wrapped_u64(rng.gen::<u64>());
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 16"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });

    let mut v = vec![F::ZERO; 1 << 20];
    for i in 0..(1 << 20) {
    	v[i] = F::from_wrapped_u64(rng.gen::<u64>());
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 20"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });
}

fn run_fft_complex<F: Field + TwoAdicField>(c: &mut Criterion)
where
    Standard: Distribution<F>,
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<F>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    let mut rng = thread_rng();

    let mut v = vec![F::ZERO; 1 << 14];
    for i in 0..(1 << 14) {
    	v[i] = rng.gen();
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 14"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });

    let mut v = vec![F::ZERO; 1 << 16];
    for i in 0..(1 << 16) {
    	v[i] = rng.gen();
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 16"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });

    let mut v = vec![F::ZERO; 1 << 20];
    for i in 0..(1 << 20) {
    	v[i] = rng.gen();
    }
    let poly = PolynomialCoeffs::new(v);
    group.bench_with_input(BenchmarkId::from_parameter("fft 1 << 20"), &poly, |bench, poly| {
        bench.iter(|| {
            poly.clone().fft()
        });
    });
}

fn run_real_fft_complex(c: &mut Criterion)
{
    let mut group = c.benchmark_group(&format!("encode::<{}>", type_name::<Mersenne31>()));
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));
    let mut rng = thread_rng();

    for size in [14, 16] {
	    let mut v1 = vec![Mersenne31::ZERO; 1 << size];
	    for i in 0..(1 << size) {
	    	v1[i] = rng.gen();
	    }
	    let mut v2 = vec![Mersenne31::ZERO; 1 << size];
	    for i in 0..(1 << size) {
	    	v2[i] = rng.gen();
	    }
	    let poly1 = PolynomialCoeffs::new(v1);
	    let poly2 = PolynomialCoeffs::new(v2);
	    group.bench_with_input(BenchmarkId::from_parameter(format!("fft 1 << {:?}", size)), &(poly1, poly2), |bench, (poly1, poly2)| {
	        bench.iter(|| {
		    	let poly = PolynomialCoeffs::<Mersenne31Complex<Mersenne31>>::new(
		    		poly1.coeffs
		    			 .clone().into_iter()
		    			 .zip(poly2.coeffs.clone().into_iter())
		    			 .map(|(a, b)| Mersenne31Complex::<Mersenne31>::new(a, b))
		    			 .collect()
		    	);
		    	let n = poly.len();
	            let values = poly.fft();
	            let values1 = (0..(n/2)).map(|k| {
	            	Mersenne31Complex::<Mersenne31>::new(
	            		(values.values[k].real() + values.values[(n-k)%n].real()).div_2exp_u64(1),
	            		(values.values[k].imag() - values.values[(n-k)%n].imag()).div_2exp_u64(1),
	            	)
	            }).collect::<Vec<Mersenne31Complex<Mersenne31>>>();
	            let values2 = (0..(n/2)).map(|k| {
	            	Mersenne31Complex::<Mersenne31>::new(
	            		(values.values[(n-k)%n].imag() + values.values[k].imag()).div_2exp_u64(1),
	            		(values.values[(n-k)%n].real() - values.values[k].real()).div_2exp_u64(1),
	            	)
	            }).collect::<Vec<Mersenne31Complex<Mersenne31>>>();
	            (values1, values2)
	        });
	    });
    }
}

criterion_group!(benches, bench_fft);
criterion_main!(benches);