use std::iter::zip;

use crate::adhesion_index::*;

#[derive(Debug, PartialEq)]
pub struct Particle {
    pub id: ParId,
    pub pos: ParPos,
    pub typ: ParticleType,
}

#[derive(serde::Deserialize, Debug)]
pub struct Particles {
    pub par_ids: Vec<ParId>,
    pub positions: Vec<f64>,
    pub types: Vec<i32>,
}

impl Particles {
    /// Find particle on index.
    pub fn at(&self, index: ParId) -> Option<Particle> {
        for (k, (id, typ)) in zip(self.par_ids.iter(), self.types.iter()).enumerate() {
            if *id == index {
                let pos = self.get_pos(k);
                return Some(Particle {
                    id: id.clone(),
                    pos: pos.clone(),
                    typ: match typ {
                        0 => ParticleType::Free,
                        1 => ParticleType::Excluded,
                        2 => ParticleType::Adhesion,
                        3 => ParticleType::Boundary,
                        _ => panic!("Unkown particle type"),
                    },
                });
            }
        }
        None
    }

    // Get positions to a ParPos, positions is stored as a 3*N array.
    fn get_pos(&self, index: usize) -> ParPos {
        (
            self.positions[3 * index],
            self.positions[3 * index + 1],
            self.positions[3 * index + 2],
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub bond_id: BondId,
    pub group: (ParId, ParId),
    pub typ: BondType,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Angle {
    pub angle_id: AngleId,
    pub group: (ParId, ParId, ParId),
    pub typ: AngleType,
}

#[derive(serde::Deserialize, Debug)]
pub struct Bonds {
    pub bond_ids: Vec<BondId>,
    pub group: Vec<(ParId, ParId)>,
    pub types: Vec<i32>,
}

impl Bonds {
    pub fn at(&self, bonds_types: &BondsTypes, index: BondId) -> Bond {
        for (k, id) in self.bond_ids.iter().enumerate() {
            if *id == index {
                let grp = self.group[k];
                let typ = self.types[k];
                return Bond {
                    bond_id: id.clone(),
                    group: grp,
                    typ: bonds_types.get(typ),
                };
            }
        }
        panic!("Particle not in data.")
    }
}

#[derive(serde::Deserialize, Debug)]
pub struct BondsTypes {
    pub bond_type_ids: Vec<i32>,
    pub k: Vec<f64>,
    pub r0: Vec<f64>,
}

impl BondsTypes {
    fn get(&self, bondid: BondId) -> BondType {
        for (id, (k, r0)) in zip(
            self.bond_type_ids.iter(),
            zip(self.k.clone(), self.r0.clone()),
        ) {
            if *id == bondid {
                return BondType::new(k, r0);
            }
        }
        panic!("Bonds Always have a type!")
    }
}

#[derive(serde::Deserialize, Debug)]
pub struct Angles {
    pub angle_cst_ids: Vec<i32>,
    pub group: Vec<(ParId, ParId, ParId)>,
    pub types: Vec<i32>,
}

impl Angles {
    pub fn at(&self, angle_types: &AnglesTypes, index: BondId) -> Angle {
        for (k, id) in self.angle_cst_ids.iter().enumerate() {
            if *id == index {
                let grp = self.group[k];
                let typ = self.types[k];
                return Angle {
                    angle_id: id.clone(),
                    group: grp,
                    typ: angle_types.get(typ),
                };
            }
        }
        panic!("Particle not in data.")
    }
}

#[derive(serde::Deserialize, Debug)]
pub struct AnglesTypes {
    pub angle_type_ids: Vec<i32>,
    pub k: Vec<f64>,
    pub t0: Vec<f64>,
}

impl AnglesTypes {
    pub fn get(&self, angle_id: AngleId) -> AngleType {
        for (id, (k, t0)) in zip(
            self.angle_type_ids.iter(),
            zip(self.k.clone(), self.t0.clone()),
        ) {
            if *id == angle_id {
                return AngleType::new(k, t0);
            }
        }
        panic!("Angles Always have a type!")
    }
}

#[derive(serde::Deserialize, Debug)]
pub struct ECMBoundary {
    pub particles: Particles,
    pub bonds: Bonds,
    pub bonds_types: BondsTypes,
    pub angles: Angles,
    pub angles_types: AnglesTypes,
}

#[cfg(test)]
mod test {
    use super::*;
    use serde_json::json;

    #[test]
    fn test() {
        let boundary = serde_json::json!({
            "particles": {
                "par_ids": [20, 50, 10, 5],
                "positions": [
                    1.0, 2.0, 3.0,
                    4.0, 5.0, 6.0,
                    7.0, 8.0, 9.0,
                    10.0, 11.0, 12.0,
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
        let boundary: ECMBoundary = serde_json::from_value(boundary).expect("Hardcoded this");
        let expected_particle = Particle {
            id: 50,
            pos: (4.0, 5.0, 6.0),
            typ: ParticleType::Free,
        };
        let particle = boundary.particles.at(50).unwrap();
        assert_eq!(particle, expected_particle);

        let expected_bond = Bond {
            bond_id: 10,
            group: (5, 20),
            typ: BondType { k: 100.0, r0: 6.25 },
        };
        let bond = boundary.bonds.at(&boundary.bonds_types, 10);
        assert_eq!(bond, expected_bond);

        let angle = boundary.angles.at(&boundary.angles_types, 15);
        let expected_angle = Angle {
            angle_id: 15,
            group: (20, 50, 10),
            typ: AngleType::new(100.0, 3.14),
        };
        assert_eq!(angle, expected_angle);
    }
}
