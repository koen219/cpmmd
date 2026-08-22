use crate::Graph;
use itertools::iproduct;

pub struct Grid3D {
    pub sizex: usize,
    pub sizey: usize,
    pub sizez: usize,
    data: Vec<u8>,
}

impl Grid3D {
    pub fn new(sizex: usize, sizey: usize, sizez: usize) -> Self {
        let data = vec![0; sizex * sizey * sizez];
        Grid3D {
            sizex,
            sizey,
            sizez,
            data,
        }
    }

    #[inline]
    pub fn from_positions(&self, x: usize, y: usize, z: usize) -> crate::Node {
        // crate::Node(x + self.sizex * y + self.sizex * self.sizey * z)
        crate::Node(x * (self.sizey * self.sizez) + y * self.sizey + z)
    }

    #[inline]
    pub fn from_node(&self, node: crate::Node) -> (i32, i32, i32) {
        let index = node.0;
        let shape = (self.sizex, self.sizey, self.sizez);
        let (depth, rem) = (index / (shape.1 * shape.2), index % (shape.1 * shape.2));
        let (height, width) = (rem / shape.2, rem % shape.2);
        (depth as i32, height as i32, width as i32)
        //        (
        //            (node.0 % self.sizex) as i32,
        //            ((node.0 % (self.sizey * self.sizex)) / self.sizex) as i32,
        //            (node.0 / (self.sizex * self.sizey)) as i32,
        //        )
    }
}

const fn twentysixneighbours() -> [(i32, i32, i32); 26] {
    let mut out = [(0, 0, 0); 26];
    let mut i = 0;
    let mut dx = -1;
    while dx <= 1 {
        let mut dy = -1;
        while dy <= 1 {
            let mut dz = -1;
            while dz <= 1 {
                if !(dx == 0 && dy == 0 && dz == 0) {
                    out[i] = (dx, dy, dz);
                    i += 1;
                }
                dz += 1;
            }
            dy += 1;
        }
        dx += 1;
    }
    out
}

const fn sixneighbours() -> [(i32, i32, i32); 6] {
    [
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    ]
}

const fn eighteenneighbours() -> [(i32, i32, i32); 18] {
    let mut out = [(0, 0, 0); 18];
    let mut i = 0;
    let mut dx: i32 = -1;
    while dx <= 1 {
        let mut dy: i32 = -1;
        while dy <= 1 {
            let mut dz: i32 = -1;
            while dz <= 1 {
                if !(dx == 0 && dy == 0 && dz == 0)
                    && !(dx.abs() == 1 && dy.abs() == 1 && dz.abs() == 1)
                {
                    out[i] = (dx, dy, dz);
                    i += 1;
                }
                dz += 1;
            }
            dy += 1;
        }
        dx += 1;
    }
    out
}

impl Graph for Grid3D {
    fn size(&self) -> usize {
        return self.sizex * self.sizey * self.sizez;
    }

    fn set(&mut self, node: crate::Node, spin: crate::Spin) {
        self.data[node.0] = spin.0 as u8;
    }

    fn get(&self, node: crate::Node) -> crate::Spin {
        crate::Spin(self.data[node.0] as i32)
    }

    fn neighbours(&self, node: crate::Node) -> crate::NodeNeighbourhood {
        let sizex_i32 = self.sizex as i32;
        let sizey_i32 = self.sizey as i32;
        let sizez_i32 = self.sizez as i32;

        let mut output = crate::NodeNeighbourhood::new(1);
        let (x_node, y_node, z_node) = self.from_node(node);
        // const NBH: [(i32, i32, i32); 26] = twentysixneighbours();
        // const NBH: [(i32, i32, i32); 6] = sixneighbours();
        const NBH: [(i32, i32, i32); 18] = eighteenneighbours();

        for (dx, dy, dz) in NBH {
            let x = x_node + dx;
            let y = y_node + dy;
            let z = z_node + dz;
            if (0..sizex_i32).contains(&x)
                && (0..sizey_i32).contains(&y)
                && (0..sizez_i32).contains(&z)
            {
                let node = self.from_positions(x as usize, y as usize, z as usize);
                output.add_neighbours(node);
            }
        }
        output
    }

    fn iter_nodes(&self) -> impl Iterator<Item = crate::Node> {
        (0..self.size()).map(|i| crate::Node(i))
    }
}

mod export {
    use super::Grid3D;
    use ndarray::Array3;
    use ndarray_npy::write_npy;
    use std::error::Error;
    use std::path::Path;

    impl Grid3D {
        pub fn export(&self, path: &Path) -> Result<(), Box<dyn Error>> {
            let shape = (self.sizex, self.sizey, self.sizez);
            let array = Array3::from_shape_vec(shape, self.data.clone())?;
            write_npy(path, &array)?;
            Ok(())
        }

        pub fn copy_data(&self) -> Vec<u8> {
            self.data.clone()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{Node, Spin};

    #[test]
    fn node_and_pos() {
        let mut grid = Grid3D::new(10, 10, 10);
        let (x, y, z) = (2, 3, 4);
        let node = grid.from_positions(x, y, z);

        assert_eq!((x as i32, y as i32, z as i32), grid.from_node(node));

        let node = Node(50);
        let (x, y, z) = grid.from_node(node);
        assert_eq!(
            node,
            grid.from_positions(x as usize, y as usize, z as usize)
        );
    }

    #[test]
    fn grid_get_and_set() {
        let mut grid = Grid3D::new(10, 10, 10);
        let (x, y, z) = (5, 5, 5);
        let node = grid.from_positions(x, y, z);
        grid.set(node, Spin(1));
        assert_eq!(grid.get(node), Spin(1));
    }

    #[test]
    fn neighbourhood() {
        let mut grid = Grid3D::new(10, 10, 10);
        let (x, y, z) = (5, 5, 5);
        let node = grid.from_positions(x, y, z);
        grid.set(node, Spin(1));

        let nbhs = grid.neighbours(node);
        for node in nbhs.neighbours.clone() {
            println!("{:?}", grid.from_node(node))
        }

        assert_eq!(nbhs.number_of_neighbours, nbhs.neighbours.len());
        assert_eq!(nbhs.number_of_neighbours, 8 + 9 + 9);
    }
}
