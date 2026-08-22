include!("generated_config.rs");
use std::cmp::max;
use std::cmp::min;
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use cpm2::grid::init2d::create_blob_in_middle;
use rand::rngs::SmallRng;
use rand::Rng;
use rand::SeedableRng;

use cpm2::cells::*;
use cpm2::energy::*;
use cpm2::graph::*;
use cpm2::grid::export2d::*;
use cpm2::grid::Grid2D;
use cpm2::ModelRules;
use cpm2::*;
struct SortingRules {
    target_area: usize,
    j_matrix: Vec<f64>,
    lambda: f64,
    output_prefix: String,
    output_interval: usize,
}

impl AdhesionConstraintParameter for SortingRules {
    fn adhesion_energy(&self, current_cell: &Cell, neighbour: &Cell) -> f64 {
        let i = max(current_cell.tau.0, neighbour.tau.0);
        let j = min(current_cell.tau.0, neighbour.tau.0);

        /*  j_matrix is symmetric and laid out in the following pattern
            j_matrix = 0
                       1 2
                       3 4 5
                       6 7 8 9
                       etc.
            Hence, the element j_matrix[i,0] is at the i * (i+1)/2 position.
        */
        let index = i * (i + 1) / 2 + j;
        self.j_matrix[index]
    }
}

impl ModelRules<Grid2D> for SortingRules {
    fn compute_energy(&mut self, grid: &Grid2D, cells: &Cells, edge: &Edge) -> f64 {
        let retracting = cells.get(edge.1);
        let extending = cells.get(edge.3);
        let mut dh = 0.0;
        dh += quadratic_constraint(
            &AreaConstraint {
                lambda: self.lambda,
                target_area: self.target_area,
            },
            &retracting,
            &extending,
        );
        dh += adhesion_constraint(grid, cells, self, edge);
        dh
    }

    fn output(&self, time: usize, grid: &Grid2D, cells: &Cells) {
        if time % self.output_interval == 0 {
            let filename = format!("{:}{:08}.png", self.output_prefix, time);
            let path = Path::new(&filename);

            // let picture = artist::parallel_draw(&grid, |spin| match cells.get(spin).tau.0 {
            let picture = draw(&grid, |spin| match cells.get(spin).tau.0 {
                1 => Color(255, 0, 0),
                2 => Color(0, 255, 0),
                _ => Color(0, 0, 0),
            });
            write_picture_as_png(&picture, path);
        }
    }
}

fn generate_cell_types(_spin: Spin, seed: u64) -> CellType {
    let mut rng = SmallRng::seed_from_u64(seed);
    let celltype: usize = rng.random_range(1..3);
    CellType(celltype)
}

pub fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() <= 1 {
        panic!("No parameter file.")
    }
    let configfile = File::open(args[1].clone()).expect("Failed to open config file");

    let config: Config =
        serde_yaml::from_reader(BufReader::new(configfile)).expect("Error parsing configuration.");

    let (grid, cells) = {
        let mut grid = Grid2D::new(config.grid_sizex, config.grid_sizey);
        let mut cells = Cells::new();
        create_blob_in_middle(
            &mut grid,
            config.init_cell_size,
            config.init_cell_number_divisions,
        );
        cells.init(&grid);
        let f: Box<dyn Fn(Spin) -> CellType> =
            Box::new(move |_spin| generate_cell_types(_spin, config.seed as u64));
        cells.set_celltypes(f);
        (grid, cells)
    };

    let rules = SortingRules {
        lambda: config.target_area_lambda,
        output_prefix: config.output_prefix.clone(),
        output_interval: config.output_interval,
        j_matrix: config.j_matrix.clone(),
        target_area: config.target_area,
    };

    let model = Model::new(
        config.seed.try_into().unwrap(),
        config.temperature,
        grid,
        cells,
        rules,
    );
    simulate(model, config.mcs);
}
