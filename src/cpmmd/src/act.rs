use cpm2::graph::*;
use std::collections::HashMap;
use std::hash::Hash;

/// The ActField stores per-node "actin" values. Internally, it is a map from `Node` to `f64`.
/// A missing entry is implicitly equivalent to a zero value.
pub struct ActField {
    values: HashMap<Node, f64>,
}

impl ActField {
    /// Create a new, empty ActField.
    pub fn new() -> Self {
        ActField {
            values: HashMap::new(),
        }
    }

    /// Retrieve the actin value at `node`. If `node` is not present, returns 0.0.
    pub fn value(&self, node: Node) -> f64 {
        *self.values.get(&node).unwrap_or(&0.0)
    }

    /// Set the actin value at `node` to `value`. If `value` ≤ 0.0, the entry is removed.
    pub fn set_value(&mut self, node: Node, value: f64) {
        if value > 0.0 {
            self.values.insert(node, value);
        } else {
            self.values.remove(&node);
        }
    }

    /// Increase the actin value at `node` by `delta`. If `node` was not present, it is inserted with initial 0.0.
    pub fn increase_value(&mut self, node: Node, delta: f64) {
        let entry = self.values.entry(node).or_insert(0.0);
        *entry += delta;
    }

    /// Decrease every actin value by 1.0. Any node whose value falls to ≤ 0.0 is removed.
    pub fn decrease(&mut self) {
        // Subtract 1.0 from each stored value...
        for v in self.values.values_mut() {
            *v -= 1.0;
        }
        // Then remove any entries that have become non-positive.
        self.values.retain(|_, &mut v| v > 0.0);
    }
}

/// Compute the geometric mean of all actin values at `node` and its neighbours that share the same spin as `node`.
///
/// If any of those actin values are ≤ 0.0, returns 0.0 immediately. If none match (impossible, because `node` itself
/// always matches its own spin), returns 0.0. Otherwise returns
///     (∏ᵢ vᵢ)^(1 / count)
/// where the product runs over all neighbours (including `node`) whose spin = main_sp.
///
/// # Panics
/// Panics if `graph.get(node)` == 0 and yet some stored actin > 0.0 (i.e. “medium” has a positive actin).
pub fn geometric_mean<G>(act_field: &ActField, graph: &G, node: Node) -> f64
where
    G: Graph,
{
    let main_spin = graph.get(node);
    // First, include the node itself:
    let mut output = 1.0_f64;
    let mut count = 0_usize;

    // Check the node itself:
    let v0 = act_field.value(node);
    if main_spin.0 == 0 {
        // If “medium” and has positive actin, that’s an error in the model.
        if v0 > 0.0 {
            panic!("geometric_mean: node in medium has positive actin value!");
        }
        // Otherwise, we treat it as zero quickly:
        return 0.0;
    }
    if v0 <= 0.0 {
        return 0.0;
    }
    count += 1;
    output *= v0;

    // Now iterate all neighbours; only include those whose spin = main_spin.
    let neighbourhood = graph.neighbours(node);
    for &nbr in &neighbourhood.neighbours {
        if graph.get(nbr) == main_spin {
            let v = act_field.value(nbr);
            if v <= 0.0 {
                return 0.0;
            }
            count += 1;
            output *= v;
        }
    }

    if count == 0 {
        // This cannot happen because `node` itself always contributes.
        return 0.0;
    }
    output.powf(1.0 / (count as f64))
}

/// Compute the ΔH contribution from the Act model, to be subtracted from the total Hamiltonian.
///
/// ΔH = (λ_act / max_act) * (GM_from − GM_to),
/// where GM_x = geometric mean of actin values around `x`.
///
/// # Panics
/// Panics if `graph.get(from)` == 0 but GM_from > 0, or if `graph.get(to)` == 0 but GM_to > 0.
/// That would mean “medium” has positive actin, which is illegal.
pub fn delta_h<G>(
    act_field: &ActField,
    graph: &G,
    from: Node,
    to: Node,
    lambda_act: f64,
    max_act: f64,
) -> f64
where
    G: Graph,
{
    let gm_from = geometric_mean(act_field, graph, from);
    if graph.get(from).0 == 0 && gm_from > 0.0 {
        panic!("delta_h: source 'from' is medium but has positive actin!");
    }
    let gm_to = geometric_mean(act_field, graph, to);
    if graph.get(to).0 == 0 && gm_to > 0.0 {
        panic!("delta_h: destination 'to' is medium but has positive actin!");
    }
    (lambda_act / max_act) * (gm_from - gm_to)
}

/// Commit a move in the Act model: if the spin at `from` > 0 (i.e. not medium),
/// then set the actin value at `to` to `max_act`. Otherwise, if spin(from) == 0,
/// explicitly set actin at `to` to 0.0 (removing any previous stored value).
pub fn commit_move<G>(act_field: &mut ActField, graph: &G, from: Node, to: Node, max_act: f64)
where
    G: Graph,
{
    let spin_from = graph.get(from);
    if spin_from.0 > 0 {
        act_field.set_value(to, max_act);
    } else {
        // If from was medium, ensure 'to' has zero actin.
        act_field.set_value(to, 0.0);
    }
}

//
// ─── UNIT TESTS ─────────────────────────────────────────────────────────────────
//

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::hash::Hash;

    //
    // A small dummy Graph implementation for testing purposes.
    //
    // We will use `Node = usize` and `Spin = i32`. Each node has a spin stored
    // in a HashMap, and we store a fixed adjacency list.
    //

    pub type Node = usize;
    pub type Spin = i32;

    impl NodeNeighbourhood {
        /// Helper constructor that begins from a Vec<Node> directly.
        pub fn from_vec(v: Vec<Node>) -> NodeNeighbourhood {
            let mut nn = NodeNeighbourhood::new(v.len());
            for &n in v.iter() {
                nn.add_neighbour(n);
            }
            nn
        }
    }

    /// DummyGraph holds:
    /// - a map `spins: HashMap<Node, Spin>`,
    /// - a map `neigh: HashMap<Node, Vec<Node>>` for adjacency,
    /// - a `node_list: Vec<Node>`.
    pub struct DummyGraph {
        spins: HashMap<Node, Spin>,
        neigh: HashMap<Node, Vec<Node>>,
        node_list: Vec<Node>,
    }

    impl DummyGraph {
        pub fn new() -> Self {
            DummyGraph {
                spins: HashMap::new(),
                neigh: HashMap::new(),
                node_list: Vec::new(),
            }
        }

        pub fn add_node(&mut self, node: Node, spin: Spin) {
            self.spins.insert(node, spin);
            self.neigh.entry(node).or_default();
            if !self.node_list.contains(&node) {
                self.node_list.push(node);
            }
        }

        /// Add an undirected edge between `a` and `b`.
        pub fn add_edge(&mut self, a: Node, b: Node) {
            self.neigh.entry(a).or_default().push(b);
            self.neigh.entry(b).or_default().push(a);
        }
    }

    impl Graph for DummyGraph {
        fn size(&self) -> usize {
            self.node_list.len()
        }

        fn set(&mut self, node: Node, spin: Spin) {
            if let Some(entry) = self.spins.get_mut(&node) {
                *entry = spin;
            }
        }

        fn get(&self, node: Node) -> Spin {
            *self.spins.get(&node).unwrap_or(&0)
        }

        fn neighbours(&self, node: Node) -> NodeNeighbourhood {
            let v = self.neigh.get(&node).cloned().unwrap_or_default();
            NodeNeighbourhood::from_vec(v)
        }

        fn iter_nodes(&self) -> Box<dyn Iterator<Item = Node> + '_> {
            Box::new(self.node_list.clone().into_iter())
        }
    }

    #[test]
    fn test_actfield_basic_operations() {
        let mut af = ActField::<Node>::new();

        // Initially, every node has value 0.0
        assert_eq!(af.value(0), 0.0);
        assert_eq!(af.value(42), 0.0);

        // Increase value at node 5 by 3.5
        af.increase_value(5, 3.5);
        assert_eq!(af.value(5), 3.5);

        // Increase again by 1.0 → total 4.5
        af.increase_value(5, 1.0);
        assert_eq!(af.value(5), 4.5);

        // Decrease all by 1.0 → node 5 ⇒ 3.5
        af.decrease();
        assert_eq!(af.value(5), 3.5);

        // Set node 5 to 0.0 explicitly → entry removed
        af.set_value(5, 0.0);
        assert_eq!(af.value(5), 0.0);

        // Decrease when empty: still no panic, stays empty
        af.decrease();
        assert_eq!(af.values.len(), 0);
    }

    #[test]
    fn test_geometric_mean_simple() {
        // Construct a tiny graph of 3 nodes: 0, 1, 2.
        // Let edges be: 0–1, 0–2, 1–2 (complete triangle).
        let mut g = DummyGraph::new();
        for node in 0..3 {
            // Let spin be 1 for nodes 0 and 1, spin 2 for node 2.
            let spin = if node < 2 { 1 } else { 2 };
            g.add_node(node, spin);
        }
        g.add_edge(0, 1);
        g.add_edge(0, 2);
        g.add_edge(1, 2);

        // Create an ActField and assign:
        //  node 0: 4.0
        //  node 1: 9.0
        //  node 2: 16.0
        let mut af = ActField::new();
        af.set_value(0, 4.0);
        af.set_value(1, 9.0);
        af.set_value(2, 16.0);

        // geometric_mean at node 0: neighbours of 0 are {1, 2}, but only those with spin=1 count (i.e. {0,1}):
        // product = 4.0 * 9.0 = 36.0, count = 2 → GM = sqrt(36.0) = 6.0.
        let gm0 = geometric_mean(&af, &g, 0);
        assert!((gm0 - 6.0).abs() < 1e-12);

        // geometric_mean at node 1: neighbours {0,2}, only spin=1 are {0,1} → same logic → 6.0
        let gm1 = geometric_mean(&af, &g, 1);
        assert!((gm1 - 6.0).abs() < 1e-12);

        // geometric_mean at node 2: spin=2, neighbours {0,1}, neither have spin=2 → but node 2 itself has value=16.0 → count=1 → GM=16.0
        let gm2 = geometric_mean(&af, &g, 2);
        assert!((gm2 - 16.0).abs() < 1e-12);

        // If we set node 2's act to 0.0, geometric_mean must return 0.0 immediately.
        af.set_value(2, 0.0);
        let gm2_zero = geometric_mean(&af, &g, 2);
        assert_eq!(gm2_zero, 0.0);

        // If spin(node)=0 and act > 0 (simulate invalid state), it should panic:
        g.set(2, 0);
        af.set_value(2, 5.0);
        let result = std::panic::catch_unwind(|| {
            geometric_mean(&af, &g, 2);
        });
        assert!(result.is_err());
    }

    #[test]
    fn test_delta_h_and_commit_move() {
        // Build a 2-node graph: 0–1. Node 0 has spin=1, node 1 has spin=2.
        let mut g = DummyGraph::new();
        g.add_node(0, 1);
        g.add_node(1, 2);
        g.add_edge(0, 1);

        let mut af = ActField::new();
        // Assign actin: node 0 → 4.0, node 1 → 9.0
        af.set_value(0, 4.0);
        af.set_value(1, 9.0);

        // λ_act = 2.0, max_act = 10.0
        let lambda_act = 2.0;
        let max_act = 10.0;
        // Compute ΔH(from=0 → to=1):
        //  GM_from = geometric_mean around 0: only 0 itself (spin=1), so GM_from=4.0
        //  GM_to   = geometric_mean around 1: only 1 itself (spin=2), so GM_to=9.0
        //  ΔH = (2 / 10) * (4 - 9) = 0.2 * (-5) = -1.0
        let dh = delta_h(&af, &g, Node(0), Node(1), lambda_act, max_act);
        assert!((dh + 1.0).abs() < 1e-12);

        // Now test the panic condition if "medium" has positive actin:
        // Set spin(1) → 0, keep act(1)=9.0
        g.set(1, 0);
        let dh_panic = std::panic::catch_unwind(|| {
            delta_h(&af, &g, Node(0), Node(1), lambda_act, max_act);
        });
        assert!(dh_panic.is_err());

        // Reset spin(1) to 2 and test commit_move:
        g.set(1, 2);
        // If spin(0)>0 (true), commit_move should set act(1)=max_act
        commit_move(&mut af, &g, Node(0), Node(1), max_act);
        assert_eq!(af.value(1), max_act);

        // Now set spin(0)=0, and commit_move: should set act(1)=0.0
        g.set(0, 0);
        commit_move(&mut af, &g, Node(0), Node(1), max_act);
        assert_eq!(af.value(1), 0.0);
    }
}
