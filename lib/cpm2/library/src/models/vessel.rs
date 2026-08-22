use std::cmp::max;
use std::cmp::min;
use std::path::Path;

use grid::export2d::Color;
use grid::init2d::create_blob_in_middle;
use grid::init2d::throw_in_cells;
use rand::rngs::SmallRng;
use rand::Rng;
use rand::SeedableRng;

use crate::cells::*;
use crate::energy::*;
use crate::grid::export2d::*;
use crate::grid::Grid2D;
use crate::*;

struct Chemotaxis {
    data: Vec<f64>,
    sizex: usize,
    sizey: usize,
}

impl Chemotaxis {
    fn new(sizex: usize, sizey: usize) -> Self {
        Chemotaxis {
            data: vec![0.0; sizex * sizey],
            sizex,
            sizey,
        }
    }

    fn decay_and_secrete(
        &mut self,
        decay_rate: f64,
        secration_rate: f64,
        delta_t: f64,
        grid: &Grid2D,
    ) {
        for node in grid.iter_nodes() {
            let Spin(spin) = grid.get(node);
            let mut delta = self.data[node.0] * decay_rate;
            if spin > 0 {
                delta += secration_rate;
            }
            self.data[node.0] += delta * delta_t;
        }
    }

    fn get(&self, node: Node) -> f64 {
        self.data[node.0]
    }

    fn diffuse(&mut self, diffusion_coefficient: f64, delta_t: f64) {
        for x in 0..self.sizex {
            self.data[0 * self.sizex + x] = 0.0;
            self.data[(self.sizey - 1) * self.sizex + x] = 0.0;
        }
        for y in 0..self.sizey {
            self.data[y * self.sizex + 0] = 0.0;
            self.data[y * self.sizex + self.sizex - 1] = 0.0;
        }
        for x in 1..(self.sizex - 1) {
            for y in 1..(self.sizey - 1) {
                let mut delta = 0.0;
                delta += self.data[y * self.sizex + x + 1];
                delta += self.data[y * self.sizex + x - 1];
                delta += self.data[(y + 1) * self.sizex + x];
                delta += self.data[(y - 1) * self.sizex + x];
                delta -= 4.0 * self.data[y * self.sizex + x];
                self.data[y * self.sizex + x] += diffusion_coefficient * delta * delta_t;
            }
        }
    }
}

struct VesselRules {
    target_area: usize,
    j_matrix: Vec<f64>,
    lambda: f64,
    output_prefix: String,
    output_interval: usize,
    chemotaxis: Chemotaxis,
}

impl AdhesionConstraintParameter for VesselRules {
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

impl ModelRules<Grid2D> for VesselRules {
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
        dh += 5.0 * (self.chemotaxis.get(edge.0) - self.chemotaxis.get(edge.2));
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

    fn after_mcs(&mut self, _time: usize, grid: &Grid2D, _cells: &Cells) {
        self.chemotaxis.decay_and_secrete(0.1, 1.0, 0.001, grid);
        for _ in 0..100 {
            self.chemotaxis.diffuse(1.0, 0.001);
        }
    }
}

pub struct VesselConfiguration {
    grid_sizex: usize,
    grid_sizey: usize,
    init_cell_size: usize,
    init_cell_number_divisions: usize,
    temperature: f64,
    target_area_lambda: f64,
    target_area: usize,
    j_matrix: Vec<f64>,
    output_interval: usize,
    output_prefix: String,
}

impl VesselConfiguration {
    pub fn new(
        grid_sizex: usize,
        grid_sizey: usize,
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

pub fn vessel(
    grid: Option<Grid2D>,
    cells: Option<Cells>,
    config: &VesselConfiguration,
    mcs: usize,
) {
    let (grid, cells) = match (grid, cells) {
        (Some(grid), Some(cells)) => (grid, cells),
        _ => {
            let mut grid = Grid2D::new(config.grid_sizex, config.grid_sizey);
            let mut cells = Cells::new();
            throw_in_cells(&mut grid, 100, config.target_area);
            cells.init(&grid);
            let f: fn(Spin) -> CellType = |_spin| CellType(1);
            cells.set_celltypes(f);
            (grid, cells)
        }
    };
    let rules = VesselRules {
        lambda: config.target_area_lambda,
        output_prefix: config.output_prefix.clone(),
        output_interval: config.output_interval,
        j_matrix: config.j_matrix.clone(),
        target_area: config.target_area,
        chemotaxis: Chemotaxis::new(config.grid_sizex, config.grid_sizey),
    };

    let model = Model::new(config.temperature, grid, cells, rules);
    // Model::new(config., config.temperature, grid, cells, rules);
    // model.cells
    simulate(model, mcs);
}
