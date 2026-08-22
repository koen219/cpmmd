use crate::*;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct Cell {
    pub spin: Spin,
    pub tau: CellType,
    pub area: usize,
    // pub target_area: usize,
    // fit_ellipse: FittingEllipse,
}

impl Cell {
    fn new(spin: Spin) -> Cell {
        let cell_type = match spin {
            Spin(0) => CellType(0),
            _ => CellType(1),
        };
        Cell {
            spin: spin,
            tau: cell_type,
            area: 0,
            //            target_area: 100,
            // fit_ellipse: FittingEllipse::new(),
        }
    }

    fn extending(&mut self) {
        self.area += 1;
    }
    fn retracting(&mut self) {
        self.area -= 1;
    }
}

pub struct Cells {
    data: HashMap<Spin, Cell>,
}

impl Cells {
    pub fn new() -> Cells {
        Cells {
            data: HashMap::new(),
        }
    }
    pub fn init<G: Graph>(&mut self, grid: &G) {
        for node in grid.iter_nodes() {
            let spin = grid.get(node);
            if let Some(cell) = self.data.get_mut(&spin) {
                cell.area += 1;
                continue;
            }
            let mut cell = Cell::new(spin);
            cell.area = 1;
            self.data.insert(spin, cell);
        }
    }

    pub fn set_celltypes(&mut self, f: Box<dyn Fn(Spin) -> CellType>) {
        for cell in self.data.values_mut() {
            if cell.spin == Spin(0) {
                cell.tau = CellType(0)
            } else {
                cell.tau = f(cell.spin);
            }
        }
    }

    pub fn commit_move(&mut self, edge: &Edge) {
        self.data.get_mut(&edge.1).unwrap().retracting();
        self.data.get_mut(&edge.3).unwrap().extending();
    }

    pub fn get(&self, spin: Spin) -> Cell {
        self.data.get(&spin).unwrap().clone()
    }

    pub fn get_mut(&mut self, spin: Spin) -> &mut Cell {
        self.data.get_mut(&spin).unwrap()
    }
}
