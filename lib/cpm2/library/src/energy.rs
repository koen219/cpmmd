use crate::{
    cells::{Cell, Cells},
    Edge, Graph,
};

pub trait TargetConstraint {
    fn target(&self, cell: &Cell) -> usize;
    fn value(&self, cell: &Cell) -> usize;
    fn scaler(&self, cell: &Cell) -> f64;
}

pub struct AreaConstraint {
    pub lambda: f64,
    pub target_area: usize,
}
impl TargetConstraint for AreaConstraint {
    fn target(&self, cell: &Cell) -> usize {
        return self.target_area;
    }
    fn value(&self, cell: &Cell) -> usize {
        return cell.area;
    }
    fn scaler(&self, cell: &Cell) -> f64 {
        return self.lambda;
    }
}

#[inline]
pub fn linear_constraint<T: TargetConstraint>(
    config: &T,
    retracting: &Cell,
    extending: &Cell,
) -> f64 {
    let mut dh = 0.0;
    if retracting.spin.0 > 0 {
        dh += -1.0 * config.scaler(retracting);
    }
    if extending.spin.0 > 0 {
        dh += 1.0 * config.scaler(extending);
    }
    dh
}

#[inline]
pub fn quadratic_constraint<T: TargetConstraint>(
    config: &T,
    retracting: &Cell,
    extending: &Cell,
) -> f64 {
    let retracting_value = config.value(retracting) as f64;
    let extending_value = config.value(extending) as f64;
    let mut dh: f64 = 0.0;
    //    let mut old: f64 = 0.0;
    //    let mut new: f64 = 0.0;
    if retracting.spin.0 > 0 {
        let target = config.target(retracting) as f64;
        let scaler = config.scaler(retracting);
        dh += scaler * (1.0 - 2.0 * retracting_value + 2.0 * target);
        //        old += (retracting_area - config.target_area(retracting) as f64).powi(2);
        //        new += (retracting_area - 1.0 - config.target_area(retracting) as f64).powi(2);
    }
    if extending.spin.0 > 0 {
        let target = config.target(extending) as f64;
        let scaler = config.scaler(extending);
        dh += scaler * (1.0 + 2.0 * extending_value - 2.0 * target);
        // old += (extending_area - config.target_area(extending) as f64).powi(2);
        // new += (extending_area + 1.0 - config.target_area(extending) as f64).powi(2);
    }
    dh
}

pub trait AdhesionConstraintParameter {
    fn adhesion_energy(&self, current_cell: &Cell, neighbour: &Cell) -> f64;
}

pub fn adhesion_constraint<G: Graph, T: AdhesionConstraintParameter>(
    grid: &G,
    cells: &Cells,
    config: &T,
    edge: &Edge,
) -> f64 {
    let mut old = 0.0f64;
    let mut new: f64 = 0.0f64;

    let retracting = cells.get(edge.1);
    let extending = cells.get(edge.3);
    for nbh in grid.neighbours(edge.0).neighbours {
        let nbh = cells.get(grid.get(nbh));
        if nbh.spin.0 != retracting.spin.0 {
            old += config.adhesion_energy(&retracting, &nbh);
        }
        if nbh.spin.0 != extending.spin.0 {
            new += config.adhesion_energy(&extending, &nbh);
        }
    }
    new - old
}
