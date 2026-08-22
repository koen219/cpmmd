use cpm2::{
    graph::{Edge, Graph, Node},
    grid::{self, Grid3D},
};
use rand::{rngs::SmallRng, seq::IteratorRandom};

use crate::{
    adhesion_index::{self, *},
    ecm_interactions::*,
};

pub struct AdhesionMover {
    pub adhesion_index: AdhesionIndex,
    suggested_move: Vec<((Node, Node), ParPos)>, // Option<((Node, Node), ParPos)>,
    annihilate_site: Option<Node>,
    pub move_adhesion_particles: MoveAdhesionParticles,
    pub remove_adhesion_particles: RemoveAdhesionParticles,
    adhesion_annihilation_penalty: f64,
    adhesion_overflow_number: usize,
    adhesion_overflow_penalty: f64,
}

impl AdhesionMover {
    pub fn new(
        adhesion_annihilation_penalty: f64,
        adhesion_overflow_number: usize,
        adhesion_overflow_penalty: f64,
    ) -> Self {
        AdhesionMover {
            adhesion_index: AdhesionIndex::new(),
            suggested_move: Vec::new(),
            annihilate_site: None,
            move_adhesion_particles: MoveAdhesionParticles::new(),
            remove_adhesion_particles: RemoveAdhesionParticles::new(),
            adhesion_annihilation_penalty,
            adhesion_overflow_number,
            adhesion_overflow_penalty,
        }
    }

    fn suggest_move(&mut self, from: Node, to: Node, delta: (f64, f64, f64)) {
        self.suggested_move.push(((from, to), delta));
    }

    fn no_move_found(&mut self, site: Node) {
        self.annihilate_site = Some(site);
    }

    pub fn accept_move(&mut self) {
        for ((from, to), delta) in self.suggested_move.iter() {
            // while let Some(((from, to), delta)) = self.suggested_move {
            for adh in self.adhesion_index.index.get(&from).unwrap().iter() {
                self.move_adhesion_particles.par_id.push(adh.par_id);
                self.move_adhesion_particles.new_pos.push((
                    adh.pos.0 + delta.0,
                    adh.pos.1 + delta.1,
                    adh.pos.2 + delta.2,
                ));
            }
            println!("Moving {:?} to {:?}", from, to);
            self.adhesion_index.move_pixel(*from, *to, *delta);
        }
        if let Some(node) = self.annihilate_site {
            for adh in self.adhesion_index.index.get(&node).unwrap().iter() {
                self.remove_adhesion_particles.par_id.push(adh.par_id);
            }
        }
    }
    pub fn reset(&mut self) {
        self.suggested_move.clear();
        self.annihilate_site = None;
    }

    fn compute_possible_retractions(&self, grid: &Grid3D, edge: Edge) -> Vec<(Node, f64)> {
        if let Some(adh_from) = self.adhesion_index.get_adhesions(edge.0) {
            return grid
                .neighbours(edge.0)
                .neighbours
                .iter()
                .filter(|node| grid.get(**node) == grid.get(edge.0))
                .map(|node| {
                    (
                        *node,
                        self.compute_energy_of_move(
                            grid,
                            grid.from_node(edge.0),
                            Some(adh_from.clone()),
                            grid.from_node(*node),
                            self.adhesion_index.get_adhesions(*node),
                        ),
                    )
                })
                .collect();
        }
        Vec::new()
    }

    fn compute_energy_of_move(
        &self,
        grid: &Grid3D,
        from: (i32, i32, i32),
        adh_from: Option<Vec<AdhesionWithEnv>>,
        to: (i32, i32, i32),
        adh_to: Option<Vec<AdhesionWithEnv>>,
    ) -> f64 {
        let mut energy: f64 = 0.0;
        let (x, y, z) = to;
        let (xp, yp, zp) = from;
        let delta = ((x - xp) as f64, (y - yp) as f64, (z - zp) as f64);

        if let Some(adhesions) = adh_from.clone() {
            energy += adhesions.iter().map(|adh| adh.energy(delta)).sum::<f64>();
        }

        let num_of_adh = match adh_from {
            None => 0,
            Some(adhesion) => adhesion.len(),
        } + match adh_to {
            None => 0,
            Some(adhesions) => adhesions.len(),
        };
        if num_of_adh > self.adhesion_overflow_number {
            energy += self.adhesion_overflow_penalty
                * (std::cmp::max(0, num_of_adh - self.adhesion_overflow_number) as f64);
        }
        energy
    }

    pub fn compute_energy(&mut self, rng: &mut SmallRng, grid: &Grid3D, edge: Edge) -> f64 {
        let mut energy = 0.0;
        // Retractions
        if let Some(_) = self.adhesion_index.get_adhesions(edge.0) {
            let min = self
                .compute_possible_retractions(grid, edge)
                .into_iter()
                .min_by(|&(_, a), &(_, b)| a.partial_cmp(&b).unwrap()); // unwrap because f64 can be nan.
            if let Some((new_pos, dh)) = min {
                energy += dh;
                let (x, y, z) = grid.from_node(edge.0);
                let (xp, yp, zp) = grid.from_node(new_pos);
                let delta = ((xp - x) as f64, (yp - y) as f64, (zp - z) as f64);
                self.suggest_move(edge.0, new_pos, delta);
            } else {
                self.no_move_found(edge.0);
                energy += self.adhesion_annihilation_penalty;
            }
        }

        // Extensions
        if let Some(adhesions) = self.adhesion_index.get_adhesions(edge.2) {
            let from = edge.2;
            let new_pos = edge.0.clone();
            let (xp, yp, zp) = grid.from_node(from);
            let (x, y, z) = grid.from_node(new_pos);
            let delta = ((x - xp) as f64, (y - yp) as f64, (z - zp) as f64);
            self.suggest_move(from, new_pos, delta);
            // Here we don't have to take into account the energy of overflow as there cannot be any fusion of adhesions on an extension.
            energy += adhesions.iter().map(|adh| adh.energy(delta)).sum::<f64>();
        }

        energy
    }
    //    pub fn compute_energy(&mut self, rng: &mut SmallRng, grid: &Grid3D, edge: Edge) -> f64 {
    //        let mut energy = 0.0;
    //
    //        // Retraction
    //        if let Some(adhesions) = self.adhesion_index.get_adhesions(edge.0) {
    //            let new_pos = grid
    //                .neighbours(edge.0)
    //                .neighbours
    //                .iter()
    //                .filter(|node| grid.get(**node) == grid.get(edge.0))
    //                .choose(rng)
    //                .cloned();
    //            if let Some(new_pos) = new_pos {
    //                let (xp, yp, zp) = grid.from_node(edge.0);
    //                let (x, y, z) = grid.from_node(new_pos);
    //                let delta = ((x - xp) as f64, (y - yp) as f64, (z - zp) as f64);
    //
    //                self.suggest_move(edge.0, new_pos, delta);
    //                energy += adhesions.iter().map(|adh| adh.energy(delta)).sum::<f64>();
    //            } else {
    //                self.no_move_found(edge.0);
    //                // No new place for the adhesion is found!
    //                energy += self.adhesion_annihilation_penalty;
    //            }
    //        }
    //        // Extension
    //        if let Some(adhesions) = self.adhesion_index.get_adhesions(edge.2) {
    //            let new_pos = edge.0.clone();
    //            let (xp, yp, zp) = grid.from_node(edge.2);
    //            let (x, y, z) = grid.from_node(new_pos);
    //            let delta = ((x - xp) as f64, (y - yp) as f64, (z - zp) as f64);
    //            self.suggest_move(edge.2, new_pos, delta);
    //            energy += adhesions.iter().map(|adh| adh.energy(delta)).sum::<f64>();
    //        }
    //
    //        // Adhesion overflow penalty
    //        if let Some(((from, to), _)) = self.suggested_move {
    //            let num_adhesion_from = match self.adhesion_index.get_adhesions(from) {
    //                Some(adhesions) => adhesions.len(),
    //                None => 0,
    //            };
    //            if num_adhesion_from > 0 {
    //                if let Some(adhesions_to) = self.adhesion_index.get_adhesions(to) {
    //                    let number_of_adhesions = num_adhesion_from + adhesions_to.len();
    //                    energy += self.adhesion_overflow_penalty
    //                        * (std::cmp::max(0, number_of_adhesions - self.adhesion_overflow_number)
    //                            as f64);
    //                }
    //            }
    //        }
    //
    //        energy
    //    }
}

#[cfg(test)]
mod test {
    use cpm2::graph::Graph;
    use cpm2::graph::Spin;
    use rand::SeedableRng;

    use super::*;
    use crate::ecm_boundary::ECMBoundary;

    fn make_boundary() -> ECMBoundary {
        let boundary = serde_json::json!({
            "particles": {
                "par_ids": [20, 50, 10, 5],
                "positions": [
                    1.0, 2.0, 3.0,
                    4.0, 5.0, 6.0, //    4.0, 5.0, 6.0,
                    7.0, 8.0, 9.0, //    7.0, 8.0, 9.0,
                    10.0, 11.0, 12.0, // 10.0, 11.0, 12.0,
                ],
                "types": [2, 0, 0, 2],
            },
            "bonds": {
                "bond_ids": [2,9,10,0],
                "group": [(20, 50), (20,10), (5,20), (5, 10)],
                "types": [0, 0, 0, 0],
            },
            "bonds_types": {
                "bond_type_ids": [0],
                "k": [100.0],
                "r0": [6.25],
            },
            "angles": {
                "angle_cst_ids": [1],
                "group": [(5, 10, 50)],
                "types": [0]
            },
            "angles_types": {
                "angle_type_ids": [0],
                "k": [100],
                "t0": [3.14],
            }

        });
        return serde_json::from_value(boundary).unwrap();
    }

    #[test]
    fn test() {
        let boundary = make_boundary();

        let mut mover = AdhesionMover::new(0.0, 0, 0.0);
        let mut grid = Grid3D::new(20, 20, 20);
        grid.set(grid.from_positions(1, 2, 3), Spin(1)); // Adhesion is here
        grid.set(grid.from_positions(10, 11, 12), Spin(1)); // Adhesion is here
        grid.set(grid.from_positions(10, 12, 12), Spin(1));
        let node1 = grid.from_positions(1, 2, 3);
        let node2 = grid.from_positions(10, 11, 12);
        let node3 = grid.from_positions(10, 12, 12);

        mover.adhesion_index.rebuild(&grid, &boundary);

        let energy = mover.compute_energy(
            &mut SmallRng::from_entropy(),
            &grid,
            Edge(node1, Spin(1), node2, Spin(0)),
        );
        assert!(energy > 0.0);
        mover.reset();

        let energy = mover.compute_energy(
            &mut SmallRng::from_entropy(),
            &grid,
            Edge(node2, Spin(1), node3, Spin(0)),
        );
        println!("{}", energy);
        assert_ne!(energy, 0.0);
        assert!(energy < 10000.0);

        let expected_move = (Some(((node2, node3), (0.0, 1.0, 0.0))));
        assert_eq!(expected_move, mover.suggested_move);
        mover.accept_move();

        dbg!(&mover.adhesion_index.index);
        dbg!(&mover.move_adhesion_particles);

        // LEts try an extension
        let node4 = grid.from_positions(10, 12, 13);
        let edge = Edge(node4, Spin(0), node3, Spin(1));
        let energy = mover.compute_energy(&mut SmallRng::from_entropy(), &grid, edge);

        println!("{:?}", mover.suggested_move);
        assert_ne!(energy, 0.0);
        // assert!(false);
    }

    //    #[test]
    //    fn test_energy() {
    //        let boundary = serde_json::from_value(serde_json::json!({
    //            "particles": {
    //                "par_ids": [1, 6, 10],
    //                "positions": [
    //                    1.0, 2.0, 3.0,
    //                    2.0, 2.0, 3.0,
    //                    1.0, 3.0, 3.0,
    //                ],
    //                "types": [0, 0, 2],
    //            },
    //            "bonds": {
    //                "bond_ids": [10, 50,],
    //                "group": [(1, 10), (6,10)],
    //                "types": [1, 2],
    //            },
    //            "bonds_types": {
    //                "bond_type_ids": [1, 2],
    //                "k": [100.0, 23.5],
    //                "r0": [1.0, 2_f64.sqrt() as f64],
    //            }
    //        }))
    //        .unwrap();
    //
    //        let mut mover = AdhesionMover::new();
    //        let mut grid = Grid3D::new(20, 20, 20);
    //
    //        let node1 = grid.from_positions(1, 2, 3);
    //        let node6 = grid.from_positions(2, 2, 3);
    //        let node10 = grid.from_positions(1, 3, 3);
    //        let new_pos = grid.from_positions(1, 3, 4);
    //        grid.set(node1, Spin(0));
    //        grid.set(node6, Spin(0));
    //        grid.set(node10, Spin(1)); // Adhesion is here
    //        grid.set(new_pos, Spin(1)); // Will be moved here
    //
    //        mover.adhesion_index.rebuild(&grid, &boundary);
    //
    //        mover.adhesion_index.get_adhesions(node10).unwrap().map(|adh| adh.de)
    //
    //        let edge = Edge(node1, Spin(0), node10, Spin(1));
    //        let energy = mover.compute_energy(&mut SmallRng::from_entropy(), &grid, edge);
    //        let expected_energy = (1.0 - 2.0_f64.sqrt()).powi(2) * 100.0 * 0.5
    //            + (2.0_f64.sqrt() - 3.0_f64.sqrt()).powi(2) * 23.5 * 0.5;
    //
    //        assert_eq!(energy, expected_energy);
    //    }
}
