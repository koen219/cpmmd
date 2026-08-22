include!("generated_config.rs");
// use generated_config::Config;

use std::cmp::max;
use std::cmp::min;
use std::default;
use std::env;

use cpm2::grid::init3d::create_blob_in_middle_3d;
use rand::rngs::SmallRng;
use rand::Rng;
use rand::SeedableRng;

use cpm2::cells::*;
use cpm2::energy::*;
use cpm2::graph::*;
use cpm2::grid::Grid3D;
use cpm2::*;

use crate::act;
use std::net::TcpStream;

use crate::adhesion_index::ParticleType;
use crate::adhesion_mover::AdhesionMover;
use crate::connected_constraint::check;
use crate::connected_constraint::BorderPixels;
use crate::ecm_boundary::ECMBoundary;
use crate::ecm_interactions::*;

fn compute_adh_zone(
    grid: &Grid3D,
    radius: f64,
    number_of_adhesions: usize,
) -> ChangeParticlesInArea {
    let spin = 1;

    let (sum_x, sum_y, sum_z, count) = grid
        .iter_nodes()
        .filter(|node| grid.get(*node).0 == spin)
        .map(|node| grid.from_node(node))
        .fold((0, 0, 0, 0), |(sx, sy, sz, count), (x, y, z)| {
            (sx + x, sy + y, sz + z, count + 1)
        });
    let avg_x = sum_x as f64 / count as f64;
    let avg_y = sum_y as f64 / count as f64;
    let avg_z = sum_z as f64 / count as f64;

    let positions = grid
        .iter_nodes()
        .filter(|node| grid.get(*node).0 == spin)
        .map(|node| grid.from_node(node))
        .filter(|(x, y, z)| {
            (*x as f64 - avg_x).powi(2) + (*y as f64 - avg_y).powi(2) + (*z as f64 - avg_z).powi(2)
                > radius.powi(2)
        });

    ChangeParticlesInArea {
        change_area: positions.collect(),
        from_type: ParticleType::Free,
        to_type: ParticleType::Adhesion,
        num_particles: number_of_adhesions,
    }
}

fn compute_excluded_zone(grid: &Grid3D, radius: f64) -> ChangeParticlesInArea {
    let spin = 1;

    let (sum_x, sum_y, sum_z, count) = grid
        .iter_nodes()
        .filter(|node| grid.get(*node).0 == spin)
        .map(|node| grid.from_node(node))
        .fold((0, 0, 0, 0), |(sx, sy, sz, count), (x, y, z)| {
            (sx + x, sy + y, sz + z, count + 1)
        });
    let avg_x = sum_x as f64 / count as f64;
    let avg_y = sum_y as f64 / count as f64;
    let avg_z = sum_z as f64 / count as f64;

    let positions: Vec<(i32, i32, i32)> = grid
        .iter_nodes()
        .filter(|node| grid.get(*node).0 == spin)
        .map(|node| grid.from_node(node))
        .filter(|(x, y, z)| {
            (*x as f64 - avg_x).powi(2) + (*y as f64 - avg_y).powi(2) + (*z as f64 - avg_z).powi(2)
                < radius.powi(2)
        })
        .collect();

    ChangeParticlesInArea {
        num_particles: positions.len(),
        change_area: positions,
        from_type: ParticleType::Free,
        to_type: ParticleType::Excluded,
    }
}

#[derive(Debug, Default)]
struct CellExtensions {
    data: Vec<(Node, Spin)>,
}

impl CellExtensions {
    fn add(&mut self, node: Node, spin: Spin) {
        self.data.push((node, spin));
    }

    fn validate_get_and_reset(&mut self, grid: &Grid3D) -> ChangeParticlesInArea {
        let change_area: Vec<(i32, i32, i32)> = self
            .data
            .iter()
            .filter_map(|(node, spin)| {
                if grid.get(*node) == *spin {
                    Some(grid.from_node(*node))
                } else {
                    None
                }
            })
            .collect();
        self.data.clear();
        ChangeParticlesInArea {
            num_particles: change_area.len(),
            change_area: change_area,
            from_type: ParticleType::Free,
            to_type: ParticleType::Adhesion,
        }
    }
}

struct SortingRules {
    config: Config,
    target_area: usize,
    j_matrix: Vec<f64>,
    check_connectivity: bool,
    stream: TcpStream,
    // adhesion_index: AdhesionIndex,
    adhesion_mover: AdhesionMover,
    rng: SmallRng,
    border_pixels: BorderPixels,
    number_of_accepted: usize,
    cell_extensions: CellExtensions,
    act_field: act::ActField,
}

impl TargetConstraint for SortingRules {
    fn target(&self, cell: &Cell) -> usize {
        if cell.spin.0 > 0 {
            return self.target_area;
        }
        0 as usize
    }
    fn value(&self, cell: &Cell) -> usize {
        if cell.spin.0 > 0 {
            return cell.area;
        }
        0 as usize
    }
    fn scaler(&self, cell: &Cell) -> f64 {
        self.config.target_area_lambda
    }
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

struct PerimiterConstrain<'a> {
    target_perimiter: usize,
    target_perimiter_lambda: f64,
    border_pixels: &'a BorderPixels,
}

impl<'a> PerimiterConstrain<'a> {
    fn new(target: usize, border: &'a BorderPixels, scaler: f64) -> Self {
        Self {
            target_perimiter: target,
            border_pixels: border,
            target_perimiter_lambda: scaler,
        }
    }
}

impl<'a> TargetConstraint for PerimiterConstrain<'a> {
    fn target(&self, _cell: &Cell) -> usize {
        self.target_perimiter
    }
    fn value(&self, cell: &Cell) -> usize {
        self.border_pixels.perimiter(cell.spin)
    }
    fn scaler(&self, _cell: &Cell) -> f64 {
        self.target_perimiter_lambda
    }
}

struct ConstantPressure {
    pressure: f64,
}

impl TargetConstraint for ConstantPressure {
    fn target(&self, _cell: &Cell) -> usize {
        0
    }
    fn value(&self, cell: &Cell) -> usize {
        cell.area
    }
    fn scaler(&self, _cell: &Cell) -> f64 {
        self.pressure
    }
}

impl ModelRules<Grid3D> for SortingRules {
    fn compute_energy(&mut self, grid: &Grid3D, cells: &Cells, edge: &Edge) -> f64 {
        let retract = cells.get(edge.1);
        let extend = cells.get(edge.3);
        let mut dh = 0.0;

        // Area constraint
        dh += (quadratic_constraint(self, &retract, &extend));
        dh += (adhesion_constraint(grid, cells, self, edge));

        let perimiter = PerimiterConstrain::new(
            self.config.target_perimiter as usize,
            &self.border_pixels,
            self.config.target_perimiter_lambda,
        );
        dh += (quadratic_constraint(&perimiter, &retract, &extend));

        let pressure = ConstantPressure {
            pressure: self.config.area_contraction,
        };
        dh += (linear_constraint(&pressure, &retract, &extend));

        if self.config.lambda_act > 0.0 {
            dh -= act::delta_h(
                &self.act_field,
                grid,
                edge.2,
                edge.0,
                self.config.lambda_act,
                self.config.max_act,
            );
        }
        if self.check_connectivity {
            let timer = std::time::Instant::now();
            dh += crate::connected_constraint::check(grid, &self.border_pixels, edge);
        }

        let timer = std::time::Instant::now();
        let dh_ecm = (self
            .adhesion_mover
            .compute_energy(&mut self.rng, grid, *edge));
        dh += dh_ecm;

        dh
    }
    fn commit_move(&mut self, grid: &Grid3D, edge: &Edge, accepted: bool) {
        if accepted {
            self.adhesion_mover.accept_move();
            if self.check_connectivity || self.config.target_perimiter_lambda > 0.0 {
                self.border_pixels.update(grid, edge);
            }
            self.number_of_accepted += 1;
            if edge.1 .0 == 0 {
                // if an extension, add it. Don't worry if it gets overwritten, we validate all extensions before sending them.
                self.cell_extensions.add(edge.0, edge.3);
            }
            act::commit_move(
                &mut self.act_field,
                grid,
                edge.2,
                edge.0,
                self.config.max_act,
            );
        }
        self.adhesion_mover.reset();
        // Nothing additional
    }

    fn output(&self, time: usize, grid: &Grid3D, cells: &Cells) {
        //        if time % self.output_interval == 0 {
        //            let filename = format!("{:}{:07}.npy", self.output_prefix, time);
        //            grid.export(Path::new(&filename)).unwrap();
        //        }
    }

    fn after_mcs(&mut self, time: usize, grid: &Grid3D, _cells: &Cells) {
        println!("Number of accepted moves {}", self.number_of_accepted);
        self.number_of_accepted = 0;
        let change_type_in_area = match self.config.adhesion_creation {
            true => Some(self.cell_extensions.validate_get_and_reset(grid)),
            false => None,
        };
        self.act_field.decrease();

        let msg = ECMInteractions {
            change_type_in_area: change_type_in_area,
            add_adhesion_particles: None,
            move_adhesion_particles: Some(std::mem::replace(
                &mut self.adhesion_mover.move_adhesion_particles,
                MoveAdhesionParticles::new(),
            )),
            remove_adhesion_particles: Some(std::mem::replace(
                &mut self.adhesion_mover.remove_adhesion_particles,
                RemoveAdhesionParticles::new(),
            )),
        };
        // let msg = msg.encode();
        let send_msg = serde_json::json!({
            "type": "STEP",
            "data": msg.encode(),
            "grid": grid.copy_data(),
            "time": time,
            // "data": {"grid": vec![0; 5]}
        });

        crate::sendrecieve::send_json_message(&mut self.stream, &send_msg)
            .expect("Message not send");
        let recieve_msg = crate::sendrecieve::receive_json_message(&mut self.stream)
            .expect("No message recieved after timestep.");
        let boundary: ECMBoundary = serde_json::from_value(recieve_msg).unwrap();
        self.adhesion_mover.adhesion_index.rebuild(grid, &boundary);

        // println!("recieved at {} {:?}", time, recieve_msg);
    }
}

pub fn sorting(grid: Option<Grid3D>, cells: Option<Cells>, config: &Config, mcs: usize) {
    println!("Starting sorting with parameters: {:?}", &config);
    let (grid, cells) = match (grid, cells) {
        (Some(grid), Some(cells)) => (grid, cells),
        _ => {
            let mut grid = Grid3D::new(config.grid_sizex, config.grid_sizey, config.grid_sizez);
            let mut cells = Cells::new();
            create_blob_in_middle_3d(
                &mut grid,
                config.init_cell_size,
                config.init_cell_number_divisions,
            );
            cells.init(&grid);
            let f: Box<dyn Fn(Spin) -> CellType> = Box::new(|_spin| {
                //                let mut rng = SmallRng::from_entropy();
                //                let celltype: usize = rng.gen_range(1..3);
                CellType(1)
            });
            cells.set_celltypes(f);
            (grid, cells)
        }
    };

    let connection_adress = env::var("CONNADRESS").expect("CONNADRESS not specified.");
    println!("Connection adress = {}", connection_adress);
    // let mut stream = TcpStream::connect("127.0.0.1:65432").unwrap();
    let mut stream = TcpStream::connect(connection_adress).expect("Connection refused?");

    let adh_zone = compute_adh_zone(
        &grid,
        config.init_adhesion_radius,
        config.num_init_adhesions,
    );

    let msg = ECMInteractions {
        change_type_in_area: Some(adh_zone),
        add_adhesion_particles: None,
        move_adhesion_particles: None,
        remove_adhesion_particles: None,
    };
    let exclusion_zone = compute_excluded_zone(&grid, config.init_adhesion_radius);
    let msg2 = ECMInteractions {
        change_type_in_area: Some(exclusion_zone),
        add_adhesion_particles: None,
        move_adhesion_particles: None,
        remove_adhesion_particles: None,
    };
    let init_message = serde_json::json!({
        "type": "INIT",
        "data": msg.encode(),
        "data2": msg2.encode(),
        "grid": grid.copy_data()
    });

    crate::sendrecieve::send_json_message(&mut stream, &init_message)
        .expect("Error in sending init message.");
    let msg = crate::sendrecieve::receive_json_message(&mut stream)
        .expect("Error in recieving init message.");
    println!("Recieved {:?}", msg);
    let boundary: ECMBoundary = serde_json::from_value(msg).unwrap();

    let mut adhesion_mover = AdhesionMover::new(
        config.adhesion_annihilation_penalty,
        config.adhesion_overflow_number,
        config.adhesion_overflow_penalty,
    );
    adhesion_mover.adhesion_index.rebuild(&grid, &boundary);

    let rules = SortingRules {
        config: config.clone(),
        j_matrix: config.j_matrix.clone(),
        target_area: config.target_area,
        check_connectivity: config.check_connectivity,
        stream: stream,
        adhesion_mover: adhesion_mover,
        rng: {
            if config.seed as u64 > 0 {
                SmallRng::seed_from_u64(config.seed as u64)
            } else {
                SmallRng::from_entropy()
            }
        },
        border_pixels: match config.check_connectivity || config.target_perimiter_lambda > 0.0 {
            true => BorderPixels::new(&grid),
            false => BorderPixels::new(&Grid3D::new(1, 1, 1)), // Creating borderpixels can take some time.
        },
        number_of_accepted: 0,
        cell_extensions: Default::default(),
        act_field: act::ActField::new(),
    };
    let seed = config.seed as u64;
    let model = Model::new(seed, config.temperature, grid, cells, rules);
    // Model::new(config., config.temperature, grid, cells, rules);
    // model.cells
    simulate(model, mcs);
}
