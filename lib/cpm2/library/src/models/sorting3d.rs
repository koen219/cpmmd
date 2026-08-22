use std::cmp::max;
use std::cmp::min;
use std::path::Path;

use crate::grid::init3d::create_blob_in_middle_3d;
use crate::grid::init3d::throw_in_cells_3d;
use rand::rngs::SmallRng;
use rand::Rng;
use rand::SeedableRng;

use crate::cells::*;
use crate::energy::*;
use crate::grid::Grid3D;
use crate::*;

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

impl ModelRules<Grid3D> for SortingRules {
    fn compute_energy(&mut self, grid: &Grid3D, cells: &Cells, edge: &Edge) -> f64 {
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
        dh += 1.0 * (grid.from_node(edge.0).2 as f64 - grid.from_node(edge.2).2 as f64);
        dh
    }

    fn output(&self, time: usize, grid: &Grid3D, cells: &Cells) {
        if time % self.output_interval == 0 {
            let filename = format!("{:}{:08}.npy", self.output_prefix, time);
            grid.export(Path::new(&filename)).unwrap();
        }
    }
}

pub struct SortingConfiguration {
    grid_sizex: usize,
    grid_sizey: usize,
    grid_sizez: usize,
    init_cell_size: usize,
    init_cell_number_divisions: usize,
    temperature: f64,
    target_area_lambda: f64,
    target_area: usize,
    j_matrix: Vec<f64>,
    output_interval: usize,
    output_prefix: String,
}

impl SortingConfiguration {
    pub fn new(
        grid_sizex: usize,
        grid_sizey: usize,
        grid_sizez: usize,
        init_cell_size: usize,
        init_cell_number_divisions: usize,
        temperature: f64,
        target_area_lambda: f64,
        target_area: usize,
        j_matrix: Vec<f64>,
        output_interval: usize,
        output_prefix: String,
    ) -> Self {
        Self {
            grid_sizex,
            grid_sizey,
            grid_sizez,
            init_cell_size,
            init_cell_number_divisions,
            temperature,
            target_area_lambda,
            target_area,
            j_matrix,
            output_interval,
            output_prefix,
        }
    }
}

pub fn sorting(
    grid: Option<Grid3D>,
    cells: Option<Cells>,
    config: &SortingConfiguration,
    mcs: usize,
) {
    let (grid, cells) = match (grid, cells) {
        (Some(grid), Some(cells)) => (grid, cells),
        _ => {
            let mut grid = Grid3D::new(config.grid_sizex, config.grid_sizey, config.grid_sizez);
            let mut cells = Cells::new();
            // throw_in_cells_3d(&mut grid, 200, config.target_area);
            create_blob_in_middle_3d(&mut grid, config.target_area, 0);
            cells.init(&grid);
            let f: fn(Spin) -> CellType = |_spin| {
                let mut rng = SmallRng::from_entropy();
                let celltype: usize = rng.gen_range(1..3);
                CellType(celltype)
            };
            cells.set_celltypes(f);
            (grid, cells)
        }
    };
    let rules = SortingRules {
        lambda: config.target_area_lambda,
        output_prefix: config.output_prefix.clone(),
        output_interval: config.output_interval,
        j_matrix: config.j_matrix.clone(),
        target_area: config.target_area,
    };

    let model = Model::new(config.temperature, grid, cells, rules);
    // Model::new(config., config.temperature, grid, cells, rules);
    // model.cells
    simulate(model, mcs);
}
