use crate::{Edge, Graph, Node};
use rand::{rngs::SmallRng, Rng};

pub trait Sampler<G: Graph> {
    fn new(grid: G) -> Self;
    fn init(&mut self);
    fn sample(&mut self, rng: &mut SmallRng) -> Option<Edge>;
    fn get_grid(&mut self) -> &mut G;
    fn into_grid(self) -> G;
    fn commit_move(&mut self, edge: &Edge) {
        self.get_grid().set(edge.0, edge.3);
    }
}

pub struct EdgeSampler<G: Graph> {
    grid: G,
    number_of_calls: usize,
}

impl<G: Graph> Sampler<G> for EdgeSampler<G> {
    fn new(grid: G) -> EdgeSampler<G> {
        EdgeSampler {
            grid,
            number_of_calls: 0,
        }
    }
    fn init(&mut self) {
        self.number_of_calls = 0;
    }

    fn sample(&mut self, rng: &mut SmallRng) -> Option<Edge> {
        loop {
            if self.number_of_calls >= self.grid.size() {
                break None;
            }
            let node = Node(rng.gen_range(0..self.grid.size()));
            let nbh = self.grid.neighbours(node);
            let nb = nbh.neighbours[rng.gen_range(0..nbh.number_of_neighbours)];
            let edge = Edge(node, self.grid.get(node), nb, self.grid.get(nb));
            self.number_of_calls += 1;
            if edge.1 == edge.3 {
                continue;
            }
            break Some(edge);
        }
    }

    fn into_grid(self) -> G {
        self.grid
    }

    fn get_grid(&mut self) -> &mut G {
        &mut self.grid
    }
}

use rand::rngs::ThreadRng;
use std::collections::HashMap;

/// EdgeSampler to handle the sampling of edges with differing spins
#[derive(Debug)]
pub struct EdgeListSampler<G: Graph> {
    grid: G,
    edgelist: HashMap<(Node, Node), usize>, // Maps each edge to a unique number
    ordered_edgelist: Vec<(Node, Node)>,    // List of edges in order for sampling
    neighbourhood_size: usize, // The size of a neighbourhood, 2D Moore -> 8 or in 3D 26
    number_of_samples: usize,
}

impl<G: Graph> Sampler<G> for EdgeListSampler<G> {
    fn new(grid: G) -> Self {
        let mut edgelist = HashMap::new();
        let mut ordered_edgelist = Vec::new();
        let mut edge_index = 0;
        let mut neighbourhood_size: usize = 1;

        for node in grid.iter_nodes() {
            let spin = grid.get(node);
            let neighbours = &grid.neighbours(node);
            neighbourhood_size = {
                let n = neighbours.number_of_neighbours;
                if n > neighbourhood_size {
                    n
                } else {
                    neighbourhood_size
                }
            };

            for &neighbour in &neighbours.neighbours {
                // Count each edge (node, neighbour) and (neighbour, node) as distinct
                if spin != grid.get(neighbour) {
                    if !edgelist.contains_key(&(node, neighbour)) {
                        edgelist.insert((node, neighbour), edge_index);
                        ordered_edgelist.push((node, neighbour));
                        edge_index += 1;
                    }

                    if !edgelist.contains_key(&(neighbour, node)) {
                        edgelist.insert((neighbour, node), edge_index);
                        ordered_edgelist.push((neighbour, node));
                        edge_index += 1;
                    }
                }
            }
        }
        dbg!(neighbourhood_size);
        dbg!(edgelist.len());
        dbg!(ordered_edgelist.len());
        Self {
            grid,
            edgelist,
            ordered_edgelist,
            neighbourhood_size,
            number_of_samples: 0,
        }
    }

    fn init(&mut self) {
        self.number_of_samples = 0;
    }

    fn sample(&mut self, rng: &mut SmallRng) -> Option<Edge> {
        while self.number_of_samples * self.neighbourhood_size < self.ordered_edgelist.len() {
            self.number_of_samples += 1;
            if let Some((node1, node2)) = self.sample_random_edge(rng) {
                return Some(Edge(
                    node1,
                    self.grid.get(node1),
                    node2,
                    self.grid.get(node2),
                ));
            }
        }
        println!("Number of samples {}", self.number_of_samples);
        None
    }

    fn get_grid(&mut self) -> &mut G {
        &mut self.grid
    }
    fn into_grid(self) -> G {
        self.grid
    }
    fn commit_move(&mut self, edge: &Edge) {
        let grid = self.get_grid();
        grid.set(edge.0, edge.3);
        self.update(edge.0);
    }
}

impl<G: Graph> EdgeListSampler<G> {
    /// Creates a new EdgeSampler from a grid

    /// Samples a random edge from the differing spin interfaces
    pub fn sample_random_edge(&mut self, rng: &mut SmallRng) -> Option<(Node, Node)> {
        if self.ordered_edgelist.is_empty() {
            None
        } else {
            let index = rng.gen_range(0..self.ordered_edgelist.len());
            Some(self.ordered_edgelist[index])
        }
    }

    // pub fn update(&mut self, grid: &G, node: Node) {
    pub fn update(&mut self, node: Node) {
        let grid = &self.grid;
        for &neighbour in &grid.neighbours(node).neighbours {
            let edge = (node, neighbour);
            let reverse_edge = (neighbour, node);

            let node_spin = grid.get(node);
            let neighbour_spin = grid.get(neighbour);

            if node_spin != neighbour_spin {
                // Add edge and its reverse if not already present
                if !self.edgelist.contains_key(&edge) {
                    self.edgelist.insert(edge, self.ordered_edgelist.len());
                    self.ordered_edgelist.push(edge);
                }
                if !self.edgelist.contains_key(&reverse_edge) {
                    self.edgelist
                        .insert(reverse_edge, self.ordered_edgelist.len());
                    self.ordered_edgelist.push(reverse_edge);
                }
            } else {
                // Remove edge and its reverse if present
                if let Some(&index) = self.edgelist.get(&edge) {
                    self.ordered_edgelist.swap_remove(index);
                    self.edgelist.remove(&edge);

                    // Update indices in edgelist after swap_remove
                    if index < self.ordered_edgelist.len() {
                        let swapped_edge = self.ordered_edgelist[index];
                        self.edgelist.insert(swapped_edge, index);
                    }
                }
                if let Some(&index) = self.edgelist.get(&reverse_edge) {
                    self.ordered_edgelist.swap_remove(index);
                    self.edgelist.remove(&reverse_edge);

                    // Update indices in edgelist after swap_remove
                    if index < self.ordered_edgelist.len() {
                        let swapped_edge = self.ordered_edgelist[index];
                        self.edgelist.insert(swapped_edge, index);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use rand::SeedableRng;

    use crate::grid::Grid2D;
    use crate::Spin;

    use super::*;

    #[test]
    fn test_empty_edge_sampler() {
        let grid = Grid2D::new(10, 10);

        let mut sampler = EdgeListSampler::new(grid);
        let mut rng = SmallRng::from_entropy();
        let result = sampler.sample_random_edge(&mut rng);
        assert_eq!(result, None);
    }

    #[test]
    fn test_adding_edge_sampler() {
        let mut grid = Grid2D::new(10, 10);
        grid.set(grid.from_positions(5 as usize, 5 as usize), Spin(1));
        let sampler = EdgeListSampler::new(grid);
        assert_eq!(sampler.ordered_edgelist.len(), 2 * 8);
    }

    #[test]
    fn test_updating_edge_sampler() {
        let mut grid = Grid2D::new(4, 3);
        grid.set(grid.from_positions(1 as usize, 1 as usize), Spin(1));
        let mut sampler = EdgeListSampler::new(grid);
        assert_eq!(sampler.ordered_edgelist.len(), 2 * 8);
        let node = sampler.grid.from_positions(2 as usize, 1 as usize);
        sampler.grid.set(node, Spin(1));
        sampler.update(node);
        // dbg!(&sampler);

        let results = [
            (0, 5),
            (1, 5),
            (2, 5),
            (4, 5),
            (8, 5),
            (9, 5),
            (10, 5),
            (1, 6),
            (2, 6),
            (3, 6),
            (7, 6),
            (11, 6),
            (10, 6),
            (9, 6),
        ]
        .map(|(i, j)| (Node(i), Node(j)));

        for result in results.iter() {
            assert_ne!(sampler.edgelist.get(result), None);
            let index = sampler.edgelist.get(result).unwrap();
            assert_eq!(sampler.ordered_edgelist[*index], *result);
        }

        assert_eq!(sampler.ordered_edgelist.len(), 2 * 14);
    }
}
