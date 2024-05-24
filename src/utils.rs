use std::{f64::consts::PI, fs::{self, File}, io::Write};

use crate::Params;

pub(crate) fn mean(R: &Vec<f64>) -> f64 {
    let mut sum = 0.0;
    for i in R {
        sum += *i as f64;
    }
    sum / R.len() as f64
}

pub(crate) fn dispersion(R: &Vec<f64>) -> f64 {
    let m = mean(R);
    let mut sum_of_squares: f64 = 0.0;
    for x in R {
        let diff = *x as f64 - m;
        sum_of_squares += diff * diff;
    }
    sum_of_squares / R.len() as f64
}


pub(crate)
// fn phi(params: &Params, x: f64) -> f64 {
//     let a = (x - params.miu) / params.sigma;
//     (1.0 / (params.sigma * (2.0 * std::f64::consts::PI).sqrt())) * (-0.5 * a * a).exp()
// } 
fn phi(params: &Params, x: f64) -> f64 {
    let SQRT_2PI: f64 = (2.0 * PI).sqrt();
    let exponent = -((x - params.miu).powf(2.0)) / (2.0 * params.sigma.powf(2.0));
    let result = (exponent).exp() / (params.sigma * SQRT_2PI);
    result // вероятность попадания в i
}
//--------------------------------------------------------------------------------
// метод прямоугольнников
fn integrate_phi(params: &Params, a: f64, b: f64, steps: usize) -> f64 {
    let mut sum = 0.0;
    let h = (b - a) / steps as f64;
    for i in 0..steps {
        let x = a + (i as f64 + 0.5) * h;
        sum += phi(params, x) * h;
    }
    sum
}


pub(crate)
fn chi_square_test(samples: &[f64], k: usize, params: &Params) -> f64 {
    let min_value = *samples.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let max_value = *samples.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    let h = (max_value - min_value) / k as f64;
    let mut n_reals = vec![0; k];

    for &sample in samples {
        let bin_index = ((sample - min_value) / h).floor() as usize;
        if bin_index < k {
            n_reals[bin_index] += 1;
        }
    }

    let N = samples.len() as f64;
    let mut chi_square = 0.0;

    for i in 0..k {
        let A = min_value + i as f64 * h;
        let B = A + h;
        let pi = integrate_phi(params, A, B, 1000);
        let n_norm = pi * N;
        let n_real = n_reals[i] as f64;
        chi_square += (n_real - n_norm).powi(2) / n_norm;
    }

    chi_square
}
//-------------------------------------------------------------------------------
fn integral_square(a: f64, i: usize, h: f64, k: usize, params: &Params) -> f64{
    let x0 = i as f64 * h + a;
    let x1 = x0+h;
    let mut sum = 0.0;

    let delta = h / k as f64; //шаг внутри отрезка; точность вычисления
    for j in 0..k {
        let x = x0 + (j as f64) * delta;
        sum += delta * phi(params, x);
    }

    sum
}

pub(crate) fn Starjess(N: u128) -> usize {
    // формула Стерджесса
    ((1.0 + 1.322 * (N as f64).log10()).floor() as usize) as usize
}

pub(crate) fn counter(a: f64, b: f64, list: &Vec<f64>) -> usize {
    list.iter().filter(|&&x| a <= x && x <= b).count()
}

pub(crate) fn full_square(A: f64, B: f64, k: usize, params: &Params) -> f64{
    let mut S = 0.0;
    let  h = (B - A) / k as f64;
    for i in 0..k {
        S += integral_square(A, i, h, k, params)
    }
    S
}

// pub(crate) fn get_real_hi2(samples: &Vec<f64>, params: &Params) -> f64 {
//     let mut hi2 = 0.0;
//     let k = Starjess(samples.len().try_into().unwrap());
//     let A = params.miu - 3.0 * params.sigma;
//     let B = params.miu + 3.0 * params.sigma;
//     let  h = (B - A) / k as f64;
//     let S_all = full_square(A, B, k, params);
//     for i in 0..k {
//         let sj = integral_square(A, i, h, k, params);
//         let pi = sj / S_all;
//         let ni_norm = pi * samples.len() as f64;
//         let ni_real = counter(A + i as f64 * h , A +(i+1)as f64 * h, &samples);

//         hi2 += (ni_real as f64 - ni_norm).powf(2.0) / ni_norm;
//     }
//     hi2
// }