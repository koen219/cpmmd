/*use pyo3::prelude::*;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn pycpm2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    Ok(())
}
*/

// use cpm2::models::sorting::sorting as inner_sorting;
// use cpm2::models::sorting::SortingConfiguration;
use cpmmd::model::sorting as inner_sorting;
use cpmmd::model::Config;
use pyo3::prelude::*;

#[pyfunction]
fn sorting(
    configfile: &str,
    //    grid_sizex: usize,
    //    grid_sizey: usize,
    //    grid_sizez: usize,
    //    init_cell_size: usize,
    //    init_cell_number_divisions: usize,
    //    temperature: f64,
    //    target_area_lambda: f64,
    //    target_area: usize,
    //    j_matrix: Vec<f64>,
    //    num_init_adhesions: i32,
    //    init_adhesion_radius: f64,
    //    adhesion_annihilation_penalty: f64,
    //    adhesion_overflow_number: usize,
    //    adhesion_overflow_penalty: f64,
    //    check_connectivity: bool,
    //    target_perimiter: f64,
    //    target_perimiter_lambda: f64,
    //    mcs: usize,
) -> PyResult<()> {
    let config: Config = serde_yaml::from_str(configfile).expect("Error parsing configuration.");
    inner_sorting(None, None, &config, config.mcs);
    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn pycpm2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sorting, m)?)?;
    Ok(())
}
