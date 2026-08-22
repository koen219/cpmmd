pub mod cells;
pub mod energy;
pub mod graph;
pub mod grid;
// pub mod models;
// mod parameter;

use cells::Cells;
use graph::edgesampler::*;
use graph::*;
use rand::{rngs::SmallRng, SeedableRng};

pub trait ModelRules<G: Graph> {
    fn compute_energy(&mut self, grid: &G, cells: &Cells, edge: &Edge) -> f64;
    fn commit_move(&mut self, _grid: &G, _edge: &Edge, _accepted: bool) {}
    fn output(&self, _time: usize, _grid: &G, _cells: &Cells) {}
    fn after_mcs(&mut self, _time: usize, _grid: &G, _cells: &Cells) {}
}

pub struct Model<G: Graph, R: ModelRules<G>> {
    pub model_rules: R,
    temperature: f64,
    cells: Cells,
    grid: G,
    rng: SmallRng,
}

impl<G: Graph, R: ModelRules<G>> Model<G, R> {
    pub fn new(seed: u64, temperature: f64, grid: G, cells: Cells, rules: R) -> Self {
        let rng = {
            if seed > 0 {
                SmallRng::seed_from_u64(seed)
            } else {
                SmallRng::from_entropy()
            }
        };
        Model {
            model_rules: rules,
            cells: cells,
            grid: grid,
            rng,
            temperature: temperature,
        }
    }

    // fn run(&mut self, grid: Grid2D, cells: Cells, max_mcs: usize) {
    fn run(mut self, max_mcs: usize) {
        let grid = self.grid;
        let mut cells = self.cells;

        let mut edge_sampler = EdgeListSampler::new(grid);

        for time in 0..max_mcs {
            println!("Time = {}", time);
            let timer = std::time::Instant::now();
            edge_sampler.init();
            while let Some(edge) = edge_sampler.sample(&mut self.rng) {
                let timer2 = std::time::Instant::now();
                let energy =
                    self.model_rules
                        .compute_energy(&edge_sampler.get_grid(), &cells, &edge);
                // println!("Energy takes {:.2?}", timer2.elapsed());
                if accept_copy_attempt(energy, self.temperature, &mut self.rng) {
                    let timer3 = std::time::Instant::now();
                    let grid = edge_sampler.get_grid();
                    self.model_rules.commit_move(grid, &edge, true);
                    // println!("Commit move takes {:.2?}", timer3.elapsed());
                    edge_sampler.commit_move(&edge);
                    cells.commit_move(&edge);
                } else {
                    let grid = edge_sampler.get_grid();
                    self.model_rules.commit_move(&grid, &edge, false);
                }
            }
            println!("One mcs takes {:.2?}", timer.elapsed());
            // grid = edge_sampler.into_grid();

            self.model_rules
                .output(time, &edge_sampler.get_grid(), &cells);
            self.model_rules
                .after_mcs(time, edge_sampler.get_grid(), &cells);
        }
    }
}

pub fn simulate<G: Graph, R: ModelRules<G>>(model: Model<G, R>, max_mcs: usize) {
    let model = model;
    model.run(max_mcs);
}
