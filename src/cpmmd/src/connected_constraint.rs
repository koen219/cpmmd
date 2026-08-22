use cpm2::{
    graph::{Edge, Graph, Node, Spin},
    grid::Grid3D,
};
use itertools::iproduct;
use std::collections::{HashMap, HashSet, VecDeque};

/// A structure for tracking border pixels of cells in a 3D grid.
pub struct BorderPixels {
    border: HashMap<i32, HashSet<Node>>,
}

/// Constructs a new `BorderPixels` instance from a `Grid3D`.
///
/// A node is considered part of a border of a cell with a given spin if it has at least one neighbor
/// with a different spin.
///
/// # Arguments
///
/// * `grid` - A reference to the 3D grid from which to compute border nodes.
///
/// # Returns
///
/// A new `BorderPixels` containing border nodes for each spin.
impl BorderPixels {
    pub fn new(grid: &Grid3D) -> Self {
        let mut border: HashMap<i32, HashSet<Node>> = HashMap::new();
        for node in grid.iter_nodes() {
            let spin = grid.get(node);
            if spin.0 == 0 {
                continue;
            }
            if grid
                .neighbours(node)
                .neighbours
                .iter()
                .any(|x| grid.get(*x) != spin)
            {
                border.entry(spin.0).or_default().insert(node);
            }
        }

        BorderPixels { border }
    }

    /// Returns the number of border nodes associated with a given spin.
    ///
    /// # Arguments
    ///
    /// * `spin` - The spin whose perimeter (border size) should be returned.
    ///
    /// # Returns
    ///
    /// The number of nodes that form the border of the given spin. Returns 0 if no border exists.
    pub fn perimiter(&self, spin: Spin) -> usize {
        match self.border.get(&spin.0) {
            Some(border) => border.len(),
            None => 0,
        }
    }

    /// Updates the border set after a change in spin at a node.
    ///
    /// This function modifies the internal border mapping to reflect a spin change defined by the given edge.
    /// It adds or removes affected nodes from the corresponding border sets.
    ///
    /// # Arguments
    ///
    /// * `grid` - A reference to the current state of the grid.
    /// * `edge` - An edge that is updated.
    ///
    pub fn update(&mut self, grid: &Grid3D, edge: &Edge) {
        let node_to_check = edge.0;
        let current_spin = edge.1;
        let new_spin = edge.3;
        let neighbourhood = grid.neighbours(node_to_check).neighbours;

        // Remove the node from the old spin set if necessary
        if let Some(set) = self.border.get_mut(&current_spin.0) {
            set.remove(&node_to_check);
            if set.is_empty() {
                self.border.remove(&current_spin.0);
            }
        }

        for &neighbour in &neighbourhood {
            let spin = grid.get(neighbour);

            // If the neighbour is now a border, ensure it's added
            if spin != new_spin {
                if spin.0 > 0 {
                    self.border.entry(spin.0).or_default().insert(neighbour);
                }
                if new_spin.0 > 0 {
                    self.border
                        .entry(new_spin.0)
                        .or_default()
                        .insert(node_to_check);
                }
            } else {
                // If the neighbour is no longer at a border, remove it
                if !grid
                    .neighbours(neighbour)
                    .neighbours
                    .iter()
                    .any(|x| grid.get(*x) != spin)
                {
                    if let Some(set) = self.border.get_mut(&spin.0) {
                        set.remove(&neighbour);
                        if set.is_empty() {
                            self.border.remove(&spin.0);
                        }
                    }
                }
            }
        }
    }

    fn is_border_connected_if_removed(&self, grid: &Grid3D, removed: &Node, spin: Spin) -> bool {
        let border_nodes = match self.border.get(&spin.0) {
            Some(border_nodes) => border_nodes,
            None => return true,
        };

        if border_nodes.len() <= 1 {
            return true; // An empty set is trivially connected
        }

        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();

        // Start from an arbitrary node in the set that is not the removed border.
        let start_node = match border_nodes.iter().find(|&&x| x != *removed) {
            Some(&node) => node,
            None => return true, // If all nodes are removed, it's trivially connected
        };

        queue.push_back(start_node);
        visited.insert(start_node);

        while let Some(node) = queue.pop_front() {
            for neighbour in grid.neighbours(node).neighbours
            //.iter()
            // .filter(|x| grid.get(**x) == spin)
            {
                if neighbour != *removed
                    && border_nodes.contains(&neighbour)
                    && visited.insert(neighbour)
                {
                    queue.push_back(neighbour);
                }
            }
        }

        // If all border nodes are visited, they are connected
        visited.len() == border_nodes.len() - 1
    }
}

/// Check if removing node breaks the local connectedness around said node.
///
/// # Arguments:
/// * `grid` - The grid data
/// * `changing_node` - The node that is about to change spin.
/// * `new_spin` - The spin which will be put in the node
///
/// # Returns
/// `true` if the neighbourhood around the pixel is connected where connected is defined with neigbhourhood size of one degree lower.
/// `false` otherwise, i.e. the neighbourhood wont be connected.
fn check_local_connectedness(grid: &Grid3D, changing_node: Node, new_spin: Spin) -> bool {
    let spin_to_check = new_spin; // grid.get(changing_node);
    let neighbours: Vec<Node> = grid
        .neighbours(changing_node)
        .neighbours
        .into_iter()
        .filter(|node| grid.get(*node) == spin_to_check)
        .collect();

    if neighbours.len() <= 1 {
        return false;
    }

    // Preform breadth-first search starting from the first neighbour.
    // Move over each pixel that is the correct spin.
    // If the remaining grid is still locally connected, then number of
    // visited nodes is the same as the number of neighbours.
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();

    let start = neighbours[0];
    visited.insert(start);
    queue.push_back(start);

    while let Some(node) = queue.pop_front() {
        // I would like to use grid.neighbours(node) here, but I need to compute the neigbhourhood one order lower.
        // I could add this into the neighbourhood function (I should at some point)
        // For now, I just copy the neighbourhood code and decrease the number.
        let (x, y, z) = grid.from_node(node);
        let neighbourhood = iproduct!(-1..=1, -1..=1, -1..=1)
            .filter(|(dx, dy, dz)| {
                !(*dx == 0 && *dy == 0 && *dz == 0)
                //    && (i32::abs(*dx) + i32::abs(*dy) + i32::abs(*dz) <= 1)
                 && (i32::abs(*dx) <= 1 && i32::abs(*dy) <= 1 && i32::abs(*dz) <= 1)
            })
            .map(|(dx, dy, dz)| (x + dx, y + dy, z + dz))
            .filter(|(xp, yp, zp)| {
                (*xp >= 0 && *yp >= 0 && *zp >= 0)
                    && (*xp as usize) < grid.sizex
                    && (*yp as usize) < grid.sizey
                    && (*zp as usize) < grid.sizez
            })
            .map(|(xp, yp, zp)| grid.from_positions(xp as usize, yp as usize, zp as usize));

        for neighbour in neighbourhood {
            // Loop over all neighbours that are not the node we are trying to remove!
            if neighbour == changing_node {
                continue;
            }
            // Neighbour should be part of the same cluster
            if grid.get(neighbour) != spin_to_check {
                continue;
            }

            if neighbours.contains(&neighbour) && !visited.contains(&neighbour) {
                visited.insert(neighbour);
                queue.push_back(neighbour);
            }
        }
    }
    visited.len() == neighbours.len()
}

pub fn check(grid: &Grid3D, border_pixels: &BorderPixels, edge: &Edge) -> f64 {
    // Extension into medium, never breaks connectedness is always ok
    let spin = edge.1;
    if spin.0 == 0 {
        return 0.0;
    }

    if border_pixels.is_border_connected_if_removed(grid, &edge.0, spin) {
        return 0.0;
    }
    return 100000.0;
}

fn check_via_local(grid: &Grid3D, edge: &Edge) -> f64 {
    // Extension
    if edge.3 .0 > 0 {
        if !(check_local_connectedness(grid, edge.0, edge.1)
            || check_local_connectedness(grid, edge.0, edge.3))
        {
            return 0.0;
        }
    }

    // Retraction
    if check_local_connectedness(grid, edge.0, edge.3)
    //== check_local_connectedness(grid, edge.0, edge.3)
    {
        return 0.0;
    }
    return 1000000.0;

    // if !((check_local_connectedness(grid, edge.0, edge.1))
    //     && (edge.1 .0 == 0 || check_local_connectedness(grid, edge.0, edge.3)))
    // // if !((edge.1 .0 > 0 || check_local_connectedness(grid, edge.0, edge.1))
    // //     && (edge.3 .0 > 0 || check_local_connectedness(grid, edge.0, edge.3)))
    // {
    //     dh += 10000.0;
    // }
    // dh
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::connected_constraint::test::graph::Spin;
    use cpm2::*;

    fn make_grid(data: Vec<(i32, i32, i32)>) -> Grid3D {
        let (dx, dy, dz) = (50, 50, 50);
        let mut grid = Grid3D::new(100, 100, 100);
        for (x, y, z) in data {
            let node = grid.from_positions((x + dx) as usize, (y + dy) as usize, (z + dz) as usize);
            grid.set(node, Spin(1));
        }
        grid
    }

    //     #[test]
    //     fn test_connected_component_count() {
    //         let plus_sign = make_grid(vec![
    //             (1, 0, 0),
    //             (-1, 0, 0),
    //             (0, 1, 0),
    //             (0, -1, 0),
    //             (0, 0, 1),
    //             (0, 0, -1),
    //             (0, 0, 0),
    //         ]);
    //         assert!(!check_local_connectedness(
    //             &plus_sign,
    //             plus_sign.from_positions(50, 50, 50),
    //             Spin(1),
    //         ));
    //         assert!(check_local_connectedness(
    //             &plus_sign,
    //             plus_sign.from_positions(49, 50, 50),
    //             Spin(1),
    //         ));
    //
    //         let square_box = make_grid(vec![
    //             (0, 0, 0),
    //             (1, 0, 0),
    //             (0, 1, 0),
    //             (1, 1, 0),
    //             (0, 0, 1),
    //             (1, 0, 1),
    //             (0, 1, 1),
    //             (1, 1, 1),
    //         ]);
    //         assert!(check_local_connectedness(
    //             &square_box,
    //             square_box.from_positions(51, 51, 51),
    //             Spin(1),
    //         ));
    //
    //         let L_shape = make_grid(vec![(0, 0, 0), (0, 1, 0), (1, 0, 0)]);
    //         assert!(check_local_connectedness(
    //             &L_shape,
    //             L_shape.from_positions(50, 50, 50),
    //             Spin(1),
    //         ),);
    //
    //         //        let line_shape = make_grid(vec![(0, 0, 0), (-1, 0, 0), (1, 0, 0)]);
    //         //        assert!(
    //         //            check_local_connectedness(&line_shape, line_shape.from_positions(50, 50, 50), Spin(1),)
    //         //                == False
    //         //        );
    //
    //         //
    //         //        let five_clusters = make_grid(vec![
    //         //            (1, 0, 0),
    //         //            (1, 1, 0),
    //         //            (-1, 0, 0),
    //         //            (0, 1, 0),
    //         //            (0, -1, 0),
    //         //            (0, 0, 1),
    //         //            (0, 0, -1),
    //         //        ]);
    //     }
    //
    #[test]
    fn test_border() {
        let grid = make_grid(vec![(0, 0, 0), (1, 0, 0), (-1, 0, 0)]);
        let border_pixels = BorderPixels::new(&grid);
        let edge = Edge(
            grid.from_positions(50, 50, 50),
            Spin(1),
            grid.from_positions(50, 50, 51),
            Spin(0),
        );
        for pixel in border_pixels.border[&1].iter() {
            println!("{:?}", grid.from_node(*pixel));
        }
        let removed_pixel = grid.from_positions(50, 50, 50);
        let result = border_pixels.is_border_connected_if_removed(&grid, &removed_pixel, Spin(1));
        assert!(!result);

        let removed_pixel = grid.from_positions(51, 50, 50);
        let result = border_pixels.is_border_connected_if_removed(&grid, &removed_pixel, Spin(1));
        assert!(result);
    }
    #[test]
    fn test_check() {
        let grid = make_grid(vec![(0, 0, 0), (1, 0, 0), (-1, 0, 0)]);
        let border_pixels = BorderPixels::new(&grid);
        let edge = Edge(
            grid.from_positions(50, 50, 50),
            Spin(1),
            grid.from_positions(50, 50, 51),
            Spin(0),
        );
        for pixel in border_pixels.border[&1].iter() {
            println!("{:?}", grid.from_node(*pixel));
        }
        println!("check = {:?}", check(&grid, &border_pixels, &edge));

        assert!(check(&grid, &border_pixels, &edge) > 0.0);

        let grid = make_grid(vec![(0, 0, 0), (1, 0, 0), (-1, 0, 0)]);
        let border_pixels = BorderPixels::new(&grid);
        let edge = Edge(
            grid.from_positions(51, 50, 50),
            Spin(1),
            grid.from_positions(50, 50, 51),
            Spin(0),
        );
        for pixel in border_pixels.border[&1].iter() {
            println!("{:?}", grid.from_node(*pixel));
        }
        assert_eq!(check(&grid, &border_pixels, &edge), 0.0);

        let grid = make_grid(vec![
            (0, 0, 0),
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (0, 0, 1),
            (0, 0, -1),
        ]);
        let border_pixels = BorderPixels::new(&grid);
        let edge = Edge(
            grid.from_positions(50, 50, 50),
            Spin(1),
            grid.from_positions(50, 50, 51),
            Spin(0),
        );
        for pixel in border_pixels.border[&1].iter() {
            println!("{:?}", grid.from_node(*pixel));
        }
        assert!(check(&grid, &border_pixels, &edge) == 0.0);
        //        let mut grid = make_grid(
        //            iproduct!(-5..=5, -5..5, -5..5)
        //                .filter(|(x, y, z)| x * x + y * y + z * z <= 9)
        //                .collect(),
        //        );
        //        grid.set(grid.from_positions(50, 50, 54), Spin(1));
        //
        //        let copy_attempt = Edge(
        //            grid.from_positions(50, 50, 53),
        //            Spin(1),
        //            grid.from_positions(0, 0, 0),
        //            Spin(0),
        //        );
        //        let border_pixels = BorderPixels::new(&grid);
        //        assert!(check(&grid, &border_pixels, &copy_attempt) > 0.0);
        //
        //        let copy_attempt = Edge(
        //            grid.from_positions(53, 50, 50),
        //            Spin(1),
        //            grid.from_positions(0, 0, 0),
        //            Spin(0),
        //        );
        //        let border_pixels = BorderPixels::new(&grid);
        //        assert!(check(&grid, &border_pixels, &copy_attempt) == 0.0);
    }

    #[test]
    fn test_border_pixels_moore() {
        let n = 2;
        let grid = make_grid(
            iproduct!(-n..=n, -n..=n, -n..=n)
                .filter(|(x, y, z)| i32::abs(*x) <= 1 && i32::abs(*y) <= 1 && i32::abs(*z) <= 1)
                .collect(),
        );
        let mut border_pixels = BorderPixels::new(&grid);
        assert_eq!(border_pixels.border[&1].len(), 26);

        let copy_attempt = Edge(
            grid.from_positions(51, 50, 50),
            Spin(1),
            grid.from_positions(52, 50, 50),
            Spin(0),
        );
        border_pixels.update(&grid, &copy_attempt);
        assert_eq!(border_pixels.border[&1].len(), 26);

        let copy_attempt = Edge(
            grid.from_positions(51, 50, 50),
            Spin(0),
            grid.from_positions(52, 50, 50),
            Spin(1),
        );
        border_pixels.update(&grid, &copy_attempt);
        assert_eq!(border_pixels.border[&1].len(), 26);

        let copy_attempt = Edge(
            grid.from_positions(52, 50, 50),
            Spin(0),
            grid.from_positions(52, 50, 50),
            Spin(1),
        );
        border_pixels.update(&grid, &copy_attempt);
        assert_eq!(border_pixels.border[&1].len(), 27);
    }
}
