use rand::{rngs::SmallRng, Rng};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Node(pub usize);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Edge(pub Node, pub Spin, pub Node, pub Spin);
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Spin(pub i32);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellType(pub usize);

pub fn accept_copy_attempt(energy: f64, temperature: f64, rng: &mut SmallRng) -> bool {
    let probability = (-energy / temperature).exp();
    let random_number = rng.gen_range(0.0..1.0);
    random_number <= probability
}

pub struct NodeNeighbourhood {
    pub number_of_neighbours: usize,
    pub neighbours: Vec<Node>,
}

impl NodeNeighbourhood {
    pub fn new(neighbourhood_degree: usize) -> NodeNeighbourhood {
        NodeNeighbourhood {
            number_of_neighbours: 0,
            neighbours: vec![],
        }
    }
    pub fn add_neighbours(&mut self, node: Node) {
        self.number_of_neighbours += 1;
        self.neighbours.push(node);
    }
}

pub trait Graph {
    fn size(&self) -> usize;
    fn set(&mut self, node: Node, spin: Spin);
    fn get(&self, node: Node) -> Spin;
    fn neighbours(&self, node: Node) -> NodeNeighbourhood;
    fn iter_nodes(&self) -> impl Iterator<Item = Node>;
}
