use crate::adhesion_index::*;

pub struct ChangeParticlesInArea {
    pub change_area: Vec<(i32, i32, i32)>,
    pub from_type: ParticleType,
    pub to_type: ParticleType,
    pub num_particles: usize,
}
impl ChangeParticlesInArea {
    pub fn new() -> Self {
        ChangeParticlesInArea {
            change_area: Vec::new(),
            from_type: ParticleType::Free,
            to_type: ParticleType::Free,
            num_particles: 0,
        }
    }
    fn encode(self) -> serde_json::Value {
        serde_json::json!({
            "change_area": self.change_area,
            "from_type": self.from_type as u8,
            "to_type": self.to_type as u8,
            "num_particles": self.num_particles,
        })
    }
}
pub struct AddAdhesionParticles;

#[derive(Debug)]
pub struct MoveAdhesionParticles {
    pub par_id: Vec<ParId>,
    pub new_pos: Vec<ParPos>,
}

impl MoveAdhesionParticles {
    pub fn new() -> Self {
        MoveAdhesionParticles {
            par_id: Vec::new(),
            new_pos: Vec::new(),
        }
    }
    fn encode(self) -> serde_json::Value {
        serde_json::json!({
            "par_id": self.par_id,
            "new_pos": self.new_pos,
        })
    }
}
#[derive(Debug)]
pub struct RemoveAdhesionParticles {
    pub par_id: Vec<ParId>,
}
impl RemoveAdhesionParticles {
    pub fn new() -> Self {
        RemoveAdhesionParticles { par_id: Vec::new() }
    }
    fn encode(self) -> serde_json::Value {
        serde_json::json!({
            "par_id": self.par_id,
        })
    }
}

pub struct ECMInteractions {
    pub change_type_in_area: Option<ChangeParticlesInArea>,
    pub add_adhesion_particles: Option<AddAdhesionParticles>,
    pub move_adhesion_particles: Option<MoveAdhesionParticles>,
    pub remove_adhesion_particles: Option<RemoveAdhesionParticles>,
}

impl ECMInteractions {
    pub fn encode(self) -> serde_json::Value {
        serde_json::json!({
            "change_type_in_area": match self.change_type_in_area {
                None => serde_json::json!({}),
                Some(change_type_in_area) => change_type_in_area.encode()
            },
            "add_adhesion_particles": match self.add_adhesion_particles {
                None => serde_json::json!({}),
                Some(add_adhesion_particles) => todo!(),
            },
            "move_adhesion_particles": match self.move_adhesion_particles {
                None => serde_json::json!({}),
                Some(move_adhesion_particles) => move_adhesion_particles.encode(),
            },
            "remove_adhesion_particles": match self.remove_adhesion_particles {
                None => serde_json::json!({}),
                Some(remove_adhesion_particles) => remove_adhesion_particles.encode(),
            },
        })
    }
}
