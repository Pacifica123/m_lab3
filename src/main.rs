pub mod utils;

use std::{f64::consts::PI, fs::{self, File}, io::Write};

use rand::Rng;

use crate::utils::{Starjess};

struct Params {
    miu : f64, // коэффициент сдвига (ср.знач)
    sigma : f64, // коэффицент масштаба (дисперсия)
}

// СТАНДАРТНЫЙ МЕТОД
// =============================================
fn gen_white_noise(params: &Params) -> f64{
    let mut rng = rand::thread_rng();
    let mut sum_tn = 0.0;
    let mut sum_tn2 = 0.0;
    let mut sum = 0.0;

    for _ in 0..12 {
        sum += rng.gen_range(0.0..1.0);
        
    }

    sum -= 6.0;

    sum
}

fn generate_standart_random(params: Params, N: usize) -> Vec<f64> {
    let mut samples = Vec::new();

    for _ in 0..N {
        let sample = params.miu + gen_white_noise(&params) * params.sigma;
        samples.push(sample);
    }

    samples
}

// РАБОТА С ФАЙЛОМ
// =============================================
fn read_parameters_from_file(file_path: &str) -> Result<(f64, f64, i32), std::io::Error> {
    let contents = fs::read_to_string(file_path)?;

    let parts: Vec<&str> = contents.split("\n\n").collect();
    let lines: Vec<&str> = parts[0].split("-----|-------|----").collect(); 
    println!("{}", lines[1]);
    let params_line = lines[1].split('|').collect::<Vec<&str>>();
    let miu: f64 = params_line[0].trim().parse().unwrap();
    let sigma: f64 = params_line[1].trim().parse().unwrap();
    let N: i32 = params_line[2].trim().parse().unwrap();

    Ok((miu, sigma, N))
}

fn get_hi2_for_s(file_path: &str, s: usize) -> Result<f64, std::io::Error> {
    let contents = fs::read_to_string(file_path)?;

    let parts: Vec<&str> = contents.split("Таблица Пирсона для alpha = 0.05").collect();
    let table_lines: Vec<&str> = parts[1].split('\n').collect();

    let mut found_table = false;
    let mut hi2 = 0.0;

    for line in table_lines {
        if found_table {
            let values: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
            if let (Some(s_val), Some(hi2_val)) = (values.get(0), values.get(1)) {
                if s_val.trim() == s.to_string() {
                    hi2 = hi2_val.trim().parse().unwrap();
                    break;
                }
            }
        } else {
            if line.starts_with('-') {
                found_table = true;
            }
        }
    }

    Ok(hi2) 
}

fn write_to_file(file_path: &str, R: &Vec<f64>, mean: f64, despersion: f64, hi2_theor: f64, hi2_real: f64) -> Result<(), std::io::Error> {
    let mut file = File::create(file_path)?;
    let content = format!("R = {:?}\n\nmean = {}\ndispersion = {}\nhi2_theor = {}    hi2_real = {}\n", R, mean, despersion, hi2_theor, hi2_real);
    file.write_all(content.as_bytes())?;
    Ok(())
}
// =============================================

fn main() {
    // тест работы с файлом
    let (miu, sigma, N) = read_parameters_from_file("src/params.txt").unwrap();
    

    let k = Starjess(N as u128); 
    let hi2_theor = get_hi2_for_s("src/params.txt", k).unwrap();

    let R_standart = generate_standart_random(Params{miu, sigma}, N as usize);
    let mean_stadart = utils::mean(&R_standart);
    let dispersion_stadart = utils::dispersion(&R_standart);
    let hi2_standart = utils::chi_square_test(&R_standart, k, &Params{miu, sigma});

    write_to_file("src/standart.txt", &R_standart, mean_stadart, dispersion_stadart, hi2_theor, hi2_standart).unwrap();

}
