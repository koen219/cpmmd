use std::collections::{hash_map::Iter, HashMap};

use crate::grid::Grid2D;
use crate::{Graph, Spin};

#[derive(Debug)]
struct InertiaTensor {
    xx: f64,
    yy: f64,
    xy: f64,
}
impl InertiaTensor {
    fn from_moments(
        area: f64,
        sum_x: f64,
        sum_y: f64,
        sum_xy: f64,
        sum_xx: f64,
        sum_yy: f64,
    ) -> InertiaTensor {
        InertiaTensor {
            xx: sum_yy - (1.0 / area) * sum_y * sum_y,
            yy: sum_xx - (1.0 / area) * sum_x * sum_x,
            xy: -(sum_xy - (1.0 / area) * sum_y * sum_x),
        }
    }

    fn largest_eigenvalue(&self) -> f64 {
        0.5 * (self.xx + self.yy)
            + 0.5
                * ((self.xx + self.yy) * (self.xx + self.yy) - 4.0 * self.xx * self.yy
                    + 4.0 * self.xy * self.xy)
                    .sqrt()
    }

    fn smallest_eigenvalue(&self) -> f64 {
        0.5 * (self.xx + self.yy)
            - 0.5
                * ((self.xx + self.yy) * (self.xx + self.yy) - 4.0 * self.xx * self.yy
                    + 4.0 * self.xy * self.xy)
                    .sqrt()
    }
}

fn solve_symmetric_degenerate_matrix(a: f64, b: f64, c: f64) -> (f64, f64) {
    let eps = 0.00001;
    // If B is super close to one of the coordinate axis return just the
    // coordinate axis.
    if b.abs() < eps {
        if a.abs() < eps {
            return (1.0, 0.0);
        }
        return (0.0, 1.0);
    }
    let minuscoverb = -c / b;
    let scale = 1.0 / (minuscoverb * minuscoverb + 1.0).sqrt();
    (minuscoverb / scale, scale)
}

#[derive(Debug, Clone, Copy)]
pub struct FittingEllipse {
    pub area: i64,
    pub sum_x: i64,
    pub sum_y: i64,
    pub sum_xy: i64,
    pub sum_xx: i64,
    pub sum_yy: i64,
}

impl FittingEllipse {
    fn new() -> FittingEllipse {
        FittingEllipse {
            area: 0,
            sum_x: 0,
            sum_y: 0,
            sum_xx: 0,
            sum_xy: 0,
            sum_yy: 0,
        }
    }

    pub fn add_site(&mut self, x: i64, y: i64) {
        self.sum_x += x;
        self.sum_y += y;
        self.sum_xx += x * x;
        self.sum_yy += y * y;
        self.sum_xy += x * y;
        self.area += 1;
    }

    pub fn remove_site(&mut self, x: i64, y: i64) {
        self.sum_x -= x;
        self.sum_y -= y;
        self.sum_xx -= x * x;
        self.sum_yy -= y * y;
        self.sum_xy -= x * y;
        self.area -= 1;
    }

    pub fn minor_axis(&self) -> (f64, f64) {
        let inertia_tensor = InertiaTensor::from_moments(
            self.area as f64,
            self.sum_x as f64,
            self.sum_y as f64,
            self.sum_xy as f64,
            self.sum_xx as f64,
            self.sum_yy as f64,
        );
        let lambda = inertia_tensor.smallest_eigenvalue();
        return solve_symmetric_degenerate_matrix(
            inertia_tensor.xx - lambda,
            inertia_tensor.xy,
            inertia_tensor.yy - lambda,
        );
    }

    pub fn major_axis(&self) -> (f64, f64) {
        let inertia_tensor = InertiaTensor::from_moments(
            self.area as f64,
            self.sum_x as f64,
            self.sum_y as f64,
            self.sum_xy as f64,
            self.sum_xx as f64,
            self.sum_yy as f64,
        );
        // Because the inertia is the higest at the small side of the ellipse.
        let lambda = inertia_tensor.largest_eigenvalue();
        return solve_symmetric_degenerate_matrix(
            inertia_tensor.xx - lambda,
            inertia_tensor.xy,
            inertia_tensor.yy - lambda,
        );
    }

    pub fn major(&self) -> f64 {
        let inertia_tensor = InertiaTensor::from_moments(
            self.area as f64,
            self.sum_x as f64,
            self.sum_y as f64,
            self.sum_xy as f64,
            self.sum_xx as f64,
            self.sum_yy as f64,
        );
        (4.0 * inertia_tensor.largest_eigenvalue() / (1.0 * self.area as f64)).sqrt()
    }

    pub fn minor(&self) -> f64 {
        let inertia_tensor = InertiaTensor::from_moments(
            self.area as f64,
            self.sum_x as f64,
            self.sum_y as f64,
            self.sum_xy as f64,
            self.sum_xx as f64,
            self.sum_yy as f64,
        );
        (4.0 * inertia_tensor.smallest_eigenvalue() / (1.0 * self.area as f64)).sqrt()
    }
    pub fn center(&self) -> (f64, f64) {
        (
            self.sum_x as f64 / self.area as f64,
            self.sum_y as f64 / self.area as f64,
        )
    }
}

#[derive(Debug)]
pub struct CellsAsEllipse {
    map: HashMap<Spin, FittingEllipse>,
}

impl CellsAsEllipse {
    pub fn new(grid: &Grid2D) -> Self {
        let mut map: HashMap<Spin, FittingEllipse> = HashMap::new();
        for x in 0..grid.sizex {
            for y in 0..grid.sizey {
                let spin = grid.get(grid.from_positions(x, y));
                if spin.0 == 0 {
                    continue;
                }
                if let Some(ellipse) = map.get_mut(&spin) {
                    ellipse.area += 1;
                    ellipse.sum_x += x as i64;
                    ellipse.sum_y += y as i64;
                    ellipse.sum_xy += x as i64 * y as i64;
                    ellipse.sum_xx += x as i64 * x as i64;
                    ellipse.sum_yy += y as i64 * y as i64;
                } else {
                    let mut ellipse = FittingEllipse::new();
                    ellipse.area += 1;
                    ellipse.sum_x += x as i64;
                    ellipse.sum_y += y as i64;
                    ellipse.sum_xy += x as i64 * y as i64;
                    ellipse.sum_xx += x as i64 * x as i64;
                    ellipse.sum_yy += y as i64 * y as i64;
                    map.insert(spin, ellipse);
                }
            }
        }
        CellsAsEllipse { map }
    }
    pub fn get(&self, spin: Spin) -> Option<&FittingEllipse> {
        self.map.get(&spin)
    }
    pub fn get_mut(&mut self, spin: Spin) -> Option<&mut FittingEllipse> {
        self.map.get_mut(&spin)
    }

    pub fn iter(&self) -> Iter<'_, Spin, FittingEllipse> {
        self.map.iter()
    }
}

#[cfg(test)]
mod tests {

    // TODO: Import test from TST.

    #[test]
    fn test() {
        assert!(true);
    }
}
