#[derive(Debug, Clone, serde::Deserialize)]
pub struct Config {
    pub grid_sizex: usize,
    pub grid_sizey: usize,
    pub init_cell_size: usize,
    pub init_cell_number_divisions: usize,
    pub temperature: f64,
    pub target_area_lambda: f64,
    pub target_area: usize,
    pub j_matrix: Vec<f64>,
    pub output_interval: usize,
    pub output_prefix: String,
    pub mcs: usize,
    pub seed: usize,
}
