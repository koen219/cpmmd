use crate::{Graph, Node, NodeNeighbourhood, Spin};

pub struct Grid2D {
    data: Vec<u8>,
    pub sizex: usize,
    pub sizey: usize,
    neighbourhood_degree: usize,
}

impl Graph for Grid2D {
    fn size(&self) -> usize {
        self.sizex * self.sizey
    }

    fn set(&mut self, node: Node, spin: Spin) {
        self.data[node.0] = spin.0 as u8;
    }

    fn get(&self, node: Node) -> Spin {
        Spin(self.data[node.0] as i32)
    }

    fn neighbours(&self, node: Node) -> NodeNeighbourhood {
        // Here I should move the 2D aspects of NodeNeighbourhood.
        let mut nbh = NodeNeighbourhood::new(self.neighbourhood_degree);

        //sizex: usize, sizey: usize, x: i32, y: i32) {
        // Precompute size bounds as i32
        let sizex_i32 = self.sizex as i32;
        let sizey_i32 = self.sizey as i32;

        let (x, y) = ((node.0 % self.sizex) as i32, (node.0 / self.sizex) as i32);

        // Use direct iteration over a static array of tuples for offsets
        const OFFSETS: [(i32, i32); 8] = [
            (-1, 0),
            (0, 1),
            (1, 0),
            (0, -1),
            (-1, 1),
            (1, 1),
            (1, -1),
            (-1, -1),
        ];

        for &(dx, dy) in &OFFSETS {
            let nx = x + dx;
            let ny = y + dy;

            // Check bounds once and only if within limits, push the neighbour
            if nx >= 0 && nx < sizex_i32 && ny >= 0 && ny < sizey_i32 {
                // Calculate index directly without intermediate conversion
                let index = (nx as usize) + (ny as usize) * self.sizex;
                nbh.add_neighbours(Node(index));
            }
        }
        nbh
    }

    fn iter_nodes(&self) -> impl Iterator<Item = Node> {
        (0..self.size()).map(|i| Node(i))
    }
}

impl Grid2D {
    pub fn new(sizex: usize, sizey: usize) -> Grid2D {
        Grid2D {
            sizex,
            sizey,
            data: vec![0; sizex * sizey],
            neighbourhood_degree: 2,
        }
    }

    pub fn from_positions<T: Into<usize>>(&self, x: T, y: T) -> Node {
        let x: usize = x.into();
        let y: usize = y.into();
        Node(x + y * self.sizex)
    }

    pub fn copy_data(&self) -> Vec<u8> {
        self.data.clone()
    }
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn create_grid() {
        let grid = Grid2D::new(20, 10);
        let sum: i32 = grid.data.iter().map(|spin| *spin as i32).sum();
        assert_eq!(sum, 0);
        assert_eq!(grid.data.len(), 20 * 10);
    }

    #[test]
    fn from_positions() {
        let grid = Grid2D::new(20, 10);
        let mut node = grid.from_positions(5u16, 0u16);
        assert_eq!(node.0, 5);
        node = grid.from_positions(5u16, 1u16);
        assert_eq!(node.0, 25);
    }

    #[test]
    fn get_and_set() {
        let mut grid = Grid2D::new(20, 10);
        let node = grid.from_positions(10u16, 4u16);
        grid.set(node, Spin(5));
        assert_eq!(grid.get(node), Spin(5));
    }

    #[test]
    fn neighbourhood_size_middle() {
        let mut grid = Grid2D::new(20, 10);
        let node = grid.from_positions(5u16, 5u16);

        let nbh = grid.neighbours(node);
        assert_eq!(nbh.number_of_neighbours, 8);
    }

    #[test]
    fn neighbourhood_edge() {
        let grid = Grid2D::new(20, 10);
        let node = grid.from_positions(5u16, 9u16);

        let nbh = grid.neighbours(node);
        assert_eq!(nbh.number_of_neighbours, 5);
    }

    #[test]
    fn neighbourhood_iterate() {
        let grid = Grid2D::new(10, 10);
        let node = grid.from_positions(5u16, 9u16);

        let nbh = grid.neighbours(node);
        let result = [
            grid.from_positions(4u16, 9u16), // 94
            grid.from_positions(6u16, 9u16), // 96
            grid.from_positions(5u16, 8u16), // 85
            grid.from_positions(6u16, 8u16), // 86
            grid.from_positions(4u16, 8u16), // 84
        ];
        for (left, right) in zip(result.iter(), nbh.neighbours.iter()) {
            assert_eq!(left, right);
        }
    }
}
