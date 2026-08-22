pub mod threedim {
    use crate::grid::ellipse;
    use crate::grid::{Grid2D, Grid3D};
    use crate::{Graph, Node, Spin};
    use rand::{thread_rng, Rng};
    fn create_cell_at_3d(
        grid: &mut Grid3D,
        middle_x: i32,
        middle_y: i32,
        middle_z: i32,
        radius: i32,
        spin: Spin,
    ) {
        for x in -radius..radius {
            for y in -radius..radius {
                for z in -radius..radius {
                    if x * x + y * y + z * z <= radius * radius {
                        grid.set(
                            grid.from_positions(
                                (middle_x + x) as usize,
                                (middle_y + y) as usize,
                                (middle_z + z) as usize,
                            ),
                            spin,
                        );
                    }
                }
            }
        }
    }
    pub fn throw_in_cells_3d(grid: &mut Grid3D, number_of_cells: usize, max_size: usize) {
        let mut rng = thread_rng();
        for spin in 1..number_of_cells + 1 {
            let radius = (max_size as f64 / 3.1415).sqrt().floor() as i32;
            let (middle_x, middle_y, middle_z) = (
                rng.gen_range(radius + 1..grid.sizex as i32 - radius - 1),
                rng.gen_range(radius + 1..grid.sizey as i32 - radius - 1),
                rng.gen_range(radius + 1..grid.sizex as i32 - radius - 1),
            );
            create_cell_at_3d(
                grid,
                middle_x,
                middle_y,
                middle_z,
                radius,
                Spin(spin as i32),
            );
        }
    }
    pub fn create_blob_in_middle_3d(
        grid: &mut Grid3D,
        size_initial_cell: usize,
        number_of_division: usize,
    ) {
        let radius = (size_initial_cell as f64 * 3.0 / (3.1415 * 4.0)).powf(1.0 / 3.0) as i32;
        // (size_initial_cell as f64 / 3.1415).sqrt().floor() as i32;
        let (middle_x, middle_y, middle_z) = (
            (grid.sizex / 2) as i32,
            (grid.sizey / 2) as i32,
            (grid.sizez / 2) as i32,
        );
        let spin = Spin(1);
        create_cell_at_3d(grid, middle_x, middle_y, middle_z, radius, spin);
    }
}

pub mod twodim {
    use crate::grid::ellipse;
    use crate::grid::{Grid2D, Grid3D};
    use crate::{Graph, Node, Spin};
    use rand::{thread_rng, Rng};
    use std::collections::HashMap;
    fn create_cell_at(grid: &mut Grid2D, middle_x: i32, middle_y: i32, radius: i32, spin: Spin) {
        for x in -radius..radius {
            for y in -radius..radius {
                if x * x + y * y <= radius * radius {
                    grid.set(
                        grid.from_positions((middle_x + x) as usize, (middle_y + y) as usize),
                        spin,
                    );
                }
            }
        }
    }

    pub fn create_blob_in_middle(
        grid: &mut Grid2D,
        size_initial_cell: usize,
        number_of_division: usize,
    ) {
        let radius = (size_initial_cell as f64 / 3.1415).sqrt().floor() as i32;
        let (middle_x, middle_y) = ((grid.sizex / 2) as i32, (grid.sizey / 2) as i32);
        let spin = Spin(1);
        create_cell_at(grid, middle_x, middle_y, radius, spin);
        let mut max_spin = 2;
        for div in 0..number_of_division {
            let ellipses = ellipse::CellsAsEllipse::new(&grid);
            let mut into_which_spin: HashMap<Spin, Spin> = HashMap::new();
            for x in 0..grid.sizex {
                for y in 0..grid.sizey {
                    let spin = grid.get(grid.from_positions(x, y));
                    if spin.0 == 0 {
                        continue;
                    }
                    let ellipse = ellipses
                        .get(spin)
                        .expect("This ellipse is always calculated.");
                    let (mx, my) = ellipse.center();
                    let division_axis = ellipse.minor_axis();
                    if (x as f64 - mx) * division_axis.1 - (y as f64 - my) * division_axis.0 > 0.0 {
                        let node = grid.from_positions(x as usize, y as usize);
                        let into_spin = match into_which_spin.get(&spin) {
                            Some(into_spin) => into_spin.clone(),
                            None => {
                                let into_spin = Spin(max_spin);
                                max_spin += 1;
                                into_which_spin.insert(spin, into_spin);
                                into_spin
                            }
                        };
                        grid.set(node, into_spin);
                    }
                }
            }
        }
    }

    pub fn randomly_populate(grid: &mut Grid2D, number_of_pixels: usize, max_spin: i32) {
        let mut rng = thread_rng();

        for _ in 0..number_of_pixels {
            let node = Node(rng.gen_range(0..grid.size()));
            let spin = Spin(rng.gen_range(0..max_spin));
            grid.set(node, spin);
        }
    }

    pub fn throw_in_cells(grid: &mut Grid2D, number_of_cells: usize, max_size: usize) {
        let mut rng = thread_rng();
        for spin in 1..number_of_cells + 1 {
            let radius = (max_size as f64 / 3.1415).sqrt().floor() as i32;
            let (middle_x, middle_y) = (
                rng.gen_range(radius + 1..grid.sizex as i32 - radius - 1),
                rng.gen_range(radius + 1..grid.sizey as i32 - radius - 1),
            );
            create_cell_at(grid, middle_x, middle_y, radius, Spin(spin as i32));
        }
    }
}
