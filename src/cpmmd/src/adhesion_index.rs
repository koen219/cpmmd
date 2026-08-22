use std::{collections::HashMap, iter::zip};

use cpm2::{graph::Node, grid::Grid3D};

use crate::ecm_boundary::ECMBoundary;

pub type ParId = i32;
pub type ParPos = (f64, f64, f64);
pub type BondId = i32;
pub type AngleId = i32;

#[derive(Copy, Clone, serde::Deserialize, Debug, PartialEq)]
pub enum ParticleType {
    Free = 0,
    Boundary = 1,
    Adhesion = 2,
    Excluded = 3,
}

#[derive(serde::Deserialize, Debug, Clone, PartialEq)]
pub struct BondType {
    pub(crate) k: f64,
    pub r0: f64,
}
impl BondType {
    pub fn new(k: f64, r0: f64) -> Self {
        BondType { k, r0 }
    }
}

#[derive(serde::Deserialize, Debug, Clone, PartialEq)]
pub struct AngleType {
    pub k: f64,
    pub t0: f64,
}
impl AngleType {
    pub fn new(k: f64, t0: f64) -> Self {
        AngleType { k, t0 }
    }
}

#[derive(Debug, Clone)]
struct AttachedBond {
    neighbour: ParPos,
    bond: BondType,
}

impl AttachedBond {
    fn energy(&self, from: ParPos, to: ParPos) -> f64 {
        let a_old = ((from.0 - self.neighbour.0).powi(2)
            + (from.1 - self.neighbour.1).powi(2)
            + (from.2 - self.neighbour.2).powi(2))
        .sqrt();
        let a_new = ((to.0 - self.neighbour.0).powi(2)
            + (to.1 - self.neighbour.1).powi(2)
            + (to.2 - self.neighbour.2).powi(2))
        .sqrt();
        let energy_old = self.bond.k * 0.5 * (a_old - self.bond.r0).powi(2);
        let energy_new = self.bond.k * 0.5 * (a_new - self.bond.r0).powi(2);
        return energy_new - energy_old;
    }
}

/// Encode the two different ways an adhesion can be linked to the rest of the simulation.
/// Situation FAR: adhesion - bead (middle) - bead (far).
/// Situation MIDDLE: bead (left) - adhesion - bead (right).
/// They are encode like this FAR(middle, far) and MIDDLE(left, right).
/// Note that for the MIDDLE type the order between left and right does not matter.
#[derive(Debug, Clone)]
enum AngularConstraintNeighbour {
    FAR(ParPos, ParPos),
    MIDDLE(ParPos, ParPos),
}
//#[derive(Debug, Clone)]
//struct AngularConstraintType {
//    k: f64,
//    theta: f64,
//}

#[derive(Debug, Clone)]
struct AngularConstraint {
    neighbours: AngularConstraintNeighbour,
    bond: AngleType,
}

impl AngularConstraint {
    fn compute_theta(left: ParPos, middle: ParPos, right: ParPos) -> f64 {
        let q1 = (left.0 - middle.0, left.1 - middle.1, left.2 - middle.2);
        let q2 = (right.0 - middle.0, right.1 - middle.1, right.2 - middle.2);
        let size_q1 = (q1.0 * q1.0 + q1.1 * q1.1 + q1.2 * q1.2).sqrt(); // powf(2);
        let size_q2 = (q2.0 * q2.0 + q2.1 * q2.1 + q2.2 * q2.2).sqrt(); // powi(2);
        let mut cos_theta = (q1.0 * q2.0 + q1.1 * q2.1 + q1.2 * q2.2) / (size_q1 * size_q2);
        if cos_theta < -1.0 {
            cos_theta = -1.0;
        }
        if cos_theta > 1.0 {
            cos_theta = 1.0
        }
        cos_theta.acos()
    }
    fn compute_energy_angle(&self, theta: f64) -> f64 {
        self.bond.k * 0.5 * (theta - self.bond.t0).powi(2)
    }

    fn energy_particle_far(&self, from: ParPos, to: ParPos, middle: ParPos, far: ParPos) -> f64 {
        let theta_old = Self::compute_theta(from, middle, far);
        let theta_new = Self::compute_theta(to, middle, far);

        self.compute_energy_angle(theta_new) - self.compute_energy_angle(theta_old)
    }

    fn energy_particle_middle(&self, from: ParPos, to: ParPos, left: ParPos, right: ParPos) -> f64 {
        let theta_old = Self::compute_theta(left, from, right);
        let theta_new = Self::compute_theta(left, to, right);

        self.compute_energy_angle(theta_new) - self.compute_energy_angle(theta_old)
    }

    pub fn energy(&self, from: ParPos, to: ParPos) -> f64 {
        match self.neighbours {
            AngularConstraintNeighbour::FAR(middle, far) => {
                self.energy_particle_far(from, to, middle, far)
            }
            AngularConstraintNeighbour::MIDDLE(left, right) => {
                self.energy_particle_middle(from, to, left, right)
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct AdhesionWithEnv {
    pub par_id: ParId,
    pub pos: ParPos,
    neighbours: Vec<AttachedBond>,
    angle_constraints: Vec<AngularConstraint>,
}

impl AdhesionWithEnv {
    pub fn new(par_id: ParId, pos: ParPos) -> Self {
        AdhesionWithEnv {
            par_id,
            pos,
            neighbours: Vec::new(),
            angle_constraints: Vec::new(),
        }
    }

    pub fn energy(&self, delta: ParPos) -> f64 {
        let from = self.pos;
        let to = (from.0 + delta.0, from.1 + delta.1, from.2 + delta.2);
        self.neighbours
            .iter()
            .map(|bond| bond.energy(from, to))
            .sum::<f64>()
            + self
                .angle_constraints
                .iter()
                .map(|angle| angle.energy(from, to))
                .sum::<f64>()
    }
}

#[derive(Debug)]
pub struct AdhesionIndex {
    pub index: HashMap<cpm2::graph::Node, Vec<AdhesionWithEnv>>,
}

impl AdhesionIndex {
    pub fn new() -> Self {
        AdhesionIndex {
            index: HashMap::new(),
        }
    }

    fn make_bond_index(boundary: &ECMBoundary) -> HashMap<ParId, Vec<BondId>> {
        // Make the map: Adhesion -> List of bond(ids). We skip all bonds that are
        // connected to an Excluded bead.
        let mut output: HashMap<ParId, Vec<BondId>> = HashMap::new();
        for (bondid, grp) in zip(boundary.bonds.bond_ids.iter(), boundary.bonds.group.iter()) {
            let p0_adh = boundary.particles.at(grp.0).unwrap().typ == ParticleType::Adhesion;
            let p0_unfit = boundary.particles.at(grp.0).unwrap().typ == ParticleType::Excluded;
            // boundary.particles.at(grp.0).unwrap().typ == ParticleType::Adhesion
            //     || boundary.particles.at(grp.0).unwrap().typ == ParticleType::Excluded;

            let p1_adh = boundary.particles.at(grp.1).unwrap().typ == ParticleType::Adhesion;
            let p1_unfit = boundary.particles.at(grp.1).unwrap().typ == ParticleType::Excluded;

            if p0_adh && !p1_unfit {
                output.entry(grp.0).or_default().push(bondid.clone());
            }
            if p1_adh && !p0_unfit {
                output.entry(grp.1).or_default().push(bondid.clone());
            }
        }
        output
    }

    fn make_angle_index(boundary: &ECMBoundary) -> HashMap<ParId, Vec<AngleId>> {
        // Make the map: Adhesion -> List of angle(ids). We skip all angles that are
        // connected to an Excluded bead.
        let mut output: HashMap<ParId, Vec<AngleId>> = HashMap::new();
        for (angleid, grp) in zip(
            boundary.angles.angle_cst_ids.iter(),
            boundary.angles.group.iter(),
        ) {
            let p0_adh = boundary.particles.at(grp.0).unwrap().typ == ParticleType::Adhesion;
            let p0_unfit = boundary.particles.at(grp.0).unwrap().typ == ParticleType::Excluded;
            // || boundary.particles.at(grp.0).unwrap().typ == ParticleType::Adhesion;

            let p1_adh = boundary.particles.at(grp.1).unwrap().typ == ParticleType::Adhesion;
            let p1_unfit = boundary.particles.at(grp.1).unwrap().typ == ParticleType::Excluded;
            // || boundary.particles.at(grp.1).unwrap().typ == ParticleType::Adhesion;

            let p2_adh = boundary.particles.at(grp.2).unwrap().typ == ParticleType::Adhesion;
            let p2_unfit = boundary.particles.at(grp.2).unwrap().typ == ParticleType::Excluded;
            // || boundary.particles.at(grp.2).unwrap().typ == ParticleType::Adhesion;

            if p0_adh && !p1_unfit && p2_unfit {
                output.entry(grp.0).or_default().push(angleid.clone());
            }
            if p1_adh && !p0_unfit && p2_unfit {
                output.entry(grp.1).or_default().push(angleid.clone());
            }
            if p2_adh && !p0_unfit && p1_unfit {
                output.entry(grp.2).or_default().push(angleid.clone());
            }
        }
        output
    }

    pub fn rebuild(&mut self, grid: &Grid3D, boundary: &ECMBoundary) {
        self.index.clear();

        let bond_index = Self::make_bond_index(boundary);
        let angle_index = Self::make_angle_index(boundary);

        for (par_id, typ) in zip(
            boundary.particles.par_ids.iter(),
            boundary.particles.types.iter(),
        ) {
            if *typ == ParticleType::Adhesion as i32 {
                let particle = boundary.particles.at(*par_id).unwrap();
                let pos = particle.pos;
                let (x, y, z) = (
                    pos.0.floor() as usize,
                    pos.1.floor() as usize,
                    pos.2.floor() as usize,
                );
                let neighbours = match bond_index.get(&par_id) {
                    Some(neighbours) => neighbours
                        .iter()
                        .map(|bid| {
                            let bond = boundary.bonds.at(&boundary.bonds_types, *bid);
                            let typ = bond.typ;
                            let grp = bond.group;

                            let neighbour = {
                                if *par_id == grp.0 {
                                    boundary.particles.at(grp.1).unwrap().pos
                                } else {
                                    boundary.particles.at(grp.0).unwrap().pos
                                }
                            };
                            AttachedBond {
                                neighbour,
                                bond: typ,
                            }
                        })
                        .collect(),
                    None => Vec::new(), // This might happen if an adhesion is clamped between two other adhesions. Then it has no bonds so it is not saved in bond_index.
                };

                let angle_constraints: Vec<AngularConstraint> = match angle_index.get(&par_id) {
                    Some(angle_cst_id) => angle_cst_id
                        .iter()
                        .map(|aid| {
                            let angle = boundary.angles.at(&boundary.angles_types, *aid);
                            let group = angle.group;

                            let neighbours = {
                                if *par_id == group.1 {
                                    AngularConstraintNeighbour::MIDDLE(
                                        boundary.particles.at(group.0).unwrap().pos,
                                        boundary.particles.at(group.2).unwrap().pos,
                                    )
                                } else {
                                    if *par_id == group.0 {
                                        AngularConstraintNeighbour::FAR(
                                            boundary.particles.at(group.1).unwrap().pos,
                                            boundary.particles.at(group.2).unwrap().pos,
                                        )
                                    } else {
                                        // par_id == group.2
                                        AngularConstraintNeighbour::FAR(
                                            boundary.particles.at(group.1).unwrap().pos,
                                            boundary.particles.at(group.0).unwrap().pos,
                                        )
                                    }
                                }
                            };
                            AngularConstraint {
                                neighbours,
                                bond: angle.typ,
                            }
                        })
                        .collect(),
                    None => Vec::new(),
                };

                let node = grid.from_positions(x, y, z);
                self.index
                    .entry(node)
                    .or_insert_with(Vec::new)
                    .push(AdhesionWithEnv {
                        par_id: *par_id,
                        pos,
                        neighbours,
                        angle_constraints,
                    });
            }
        }
    }

    pub fn move_pixel(&mut self, from: Node, to: Node, difference: ParPos) {
        if let Some(adhesion_in_extending_pixel) = self.index.get(&from).cloned() {
            // Process and move adhesions
            self.index.entry(to).or_insert_with(Vec::new).extend(
                adhesion_in_extending_pixel.into_iter().map(|mut adhev| {
                    adhev.pos.0 += difference.0;
                    adhev.pos.1 += difference.1;
                    adhev.pos.2 += difference.2;
                    adhev // Return the modified item
                }),
            );

            // Remove the original entry
            self.index.remove(&from);
        }
    }

    fn remove_adhesion(&mut self, pixel: Node) {
        // This doesn not communicate removal of an adhesion, it just updates the index.
        self.index.remove(&pixel);
    }

    pub fn get_adhesions(&self, node: Node) -> Option<Vec<AdhesionWithEnv>> {
        // fn is_adhesion(&self, node: cpm2::graph::Node) -> bool {
        return self.index.get(&node).cloned();
    }
}

#[cfg(test)]
mod test {
    use cpm2::graph::Graph;
    use cpm2::graph::Spin;

    use crate::ecm_boundary::Bond;

    use super::*;

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
                "angle_cst_ids": [15, 2],
                "group": [(20, 50, 10), (10, 5, 50)],
                "types": [0, 0]
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

        let mut index = AdhesionIndex::new();
        let mut grid = Grid3D::new(20, 20, 20);
        grid.set(grid.from_positions(1, 2, 3), Spin(1));
        grid.set(grid.from_positions(10, 11, 12), Spin(1));
        grid.set(grid.from_positions(10, 12, 12), Spin(1));

        index.rebuild(&grid, &boundary);

        assert_eq!(index.index.len(), 2);
        let node1 = grid.from_positions(1, 2, 3);
        let node2 = grid.from_positions(10, 11, 12);
        assert!(index.index.keys().all(|k| *k == node1 || *k == node2));

        let adhe = index.index[&node1].clone();
        assert_eq!(adhe.len(), 1);
        let adhe = &adhe[0];

        assert_eq!(adhe.par_id, 20);
        assert_eq!(adhe.pos, (1.0, 2.0, 3.0));
        assert_eq!(adhe.neighbours.len(), 2);
        assert!(adhe
            .neighbours
            .iter()
            .all(|ab| ab.neighbour == (4.0, 5.0, 6.0) || ab.neighbour == (7.0, 8.0, 9.0)));

        let adhe = index.index[&node2].clone();
        assert_eq!(adhe.len(), 1);
        let adhe = adhe[0].clone();
        assert_eq!(adhe.par_id, 5);
        assert_eq!(adhe.pos, (10.0, 11.0, 12.0));
        assert_eq!(adhe.neighbours.len(), 1);
        assert_eq!(adhe.neighbours[0].neighbour, ((7.0, 8.0, 9.0)));

        //        assert_eq!(adhe.angle_constraints.len(), 1);
        //        assert_eq!(
        //            adhe.angle_constraints[0],
        //            AngularConstraint {
        //                neighbours: AngluarConstraintNeighbour {}
        //
        //            }
        //        );
    }

    #[test]
    fn energy() {
        let bond1 = AttachedBond {
            neighbour: (1.0, 2.0, 3.0),
            bond: BondType { k: 100.0, r0: 1.0 },
        };
        let bond2 = AttachedBond {
            neighbour: (2.0, 2.0, 3.0),
            bond: BondType {
                k: 23.0,
                r0: 2.0_f64.sqrt(),
            },
        };

        let mut adh = AdhesionWithEnv {
            par_id: 1,
            pos: (1.0, 3.0, 3.0),
            neighbours: vec![bond1, bond2],
            angle_constraints: Vec::new(),
        };

        assert_eq!(adh.energy((0.0, 0.0, 0.0)), 0.0);

        let old_energy = 0.0;
        let delta_energy = (2.0_f64.sqrt() - 1.0).powi(2) * 100.0 * 0.5
            + (2.0_f64.sqrt() - 3.0_f64.sqrt()).powi(2) * 23.0 * 0.5;
        assert_eq!(adh.energy((0.0, 0.0, 1.0)), delta_energy);
    }
}
