use std::collections::HashMap;
use std::hash::Hash;
use num_bigint::{BigUint, ToBigUint};
use num::Integer;
use crate::error::Error;

#[derive(Debug)]
pub enum FermiStats {
  Boson, Fermion
}
#[derive(Debug)]
pub struct Flavor {
  pub stats: FermiStats,
  pub charged: bool,
}
// degree (in,out) for one flavor at a vertex
// (note: for uncharged flavors, out=0)
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct Degree {
  pub i: usize,
  pub o: usize,
}
impl Degree {
  pub fn new() -> Self {
    Self { i: 0, o: 0 }
  }
  pub fn conj(&self) -> Self {
    Self { i: self.o, o: self.i }
  }
}
#[derive(Debug)]
pub struct VertexKind {
  pub degrees: Vec<Degree>,
}
#[derive(Debug)]
pub struct VertexState {
  occupied: Vec<Degree>,
  kind: usize,
}

#[derive(Debug)]
pub struct Theory {
  pub flavors: Vec<Flavor>,
  pub vertices: Vec<VertexKind>,
}
impl Theory {
  pub fn make_vertex(&self, kind: usize) -> VertexState {
    let n_flavor = self.flavors.len();
    assert_eq!(self.vertices[kind].degrees.len(), n_flavor);
    VertexState {
      occupied: vec![Degree::new(); n_flavor],
      kind
    }
  }
}

#[derive(Debug,Clone,PartialEq)]
pub struct Bond {
  pub flavor: usize,
  pub degree: Degree,
  pub i: usize,
  pub j: usize,
}
#[derive(Debug)]
pub struct Contraction {
  pub bonds: Vec<Bond>,
}

fn partition_into(
  n: usize, sizes: &[usize])
  -> Vec<Vec<usize>> {
  if n == 0 {
    return vec![vec![0; sizes.len()]];
  }
  if sizes.len() == 0 {
    return vec![];
  }
  let next_size = sizes[sizes.len()-1];
  let mut out = vec![];
  for m in 0..=n {
    if m > next_size {
      break;
    }
    for mut soln in partition_into(n-m, &sizes[..sizes.len()-1]) {
      soln.push(m);
      out.push(soln);
    }
  }
  out
}

pub fn enumerate_contractions(
  vertices: &mut Vec<VertexState>, i: usize, a: usize, theory: &Theory)
  -> Vec<Contraction> {
  // println!("args {:?} {} {}", vertices, i, a);
  // base case 1
  if i == vertices.len() {
    return vec![Contraction { bonds: vec![] } ];
  }

  let state = &vertices[i];
  let kind = &theory.vertices[state.kind];

  // base case 2
  if a == state.occupied.len() {
    return enumerate_contractions(vertices, i+1, 0, theory);
  }

  // recursive case
  assert!(kind.degrees[a].i >= state.occupied[a].i);
  assert!(kind.degrees[a].o >= state.occupied[a].o);
  let is_charged = theory.flavors[a].charged;
  if !is_charged {
    assert_eq!(kind.degrees[a].o, 0);
  }
  let (n_out, n_in) = if is_charged {
    // for charged wire in <- out and out -> in
    (kind.degrees[a].o - state.occupied[a].o,
     kind.degrees[a].i - state.occupied[a].i)
  }
  else {
    // for uncharged wire in <-> in and ignore out
    // (we store it anyway for uniformity)
    (kind.degrees[a].i - state.occupied[a].i,
     kind.degrees[a].o - state.occupied[a].o)
  };
  let mut n_other_in = vec![];
  let mut n_other_out = vec![];
  for vp in vertices[i+1..].iter() {
    let kind_p = &theory.vertices[vp.kind];
    assert!(kind_p.degrees[a].i >= vp.occupied[a].i);
    assert!(kind_p.degrees[a].o >= vp.occupied[a].o);
    n_other_in.push(kind_p.degrees[a].i - vp.occupied[a].i);
    n_other_out.push(kind_p.degrees[a].o - vp.occupied[a].o);
    if !is_charged {
      assert_eq!(vp.occupied[a].o, 0);
      assert_eq!(kind_p.degrees[a].o, 0);
    }
  }
  // println!("Assigning for {} {} -> {:?} {:?}", open_x, open_y, other_open_x, other_open_y);

  let assignments_out = partition_into(n_out, &n_other_in);
  let assignments_in = partition_into(n_in, &n_other_out);
  let mut out = vec![];
  for assign_out in assignments_out.iter() {
    assert_eq!(assign_out.len(), vertices.len() - i - 1);
    for j in i+1..vertices.len() {
      vertices[j].occupied[a].i += assign_out[j-i-1];
    }
    for assign_in in assignments_in.iter() {
      assert_eq!(assign_in.len(), vertices.len() - i - 1);
      for j in i+1..vertices.len() {
        vertices[j].occupied[a].o += assign_in[j-i-1];
      }
      let mut contraction_i = vec![];
      for j in i+1..vertices.len() {
        let degree = if is_charged {
          Degree { i: assign_in[j-i-1], o: assign_out[j-i-1] }
        }
        else {
          assert_eq!(assign_in[j-i-1], 0);
          Degree { i: assign_out[j-i-1], o: 0 }
        };
        contraction_i.push(Bond {flavor: a, degree, i, j});
      }
      // println!("Recursing... ");
      let sub = enumerate_contractions(vertices, i, a+1, theory);
      // println!("... result {:?}", sub);
      for mut sub_contraction in sub {
        sub_contraction.bonds.extend(contraction_i.clone());
        out.push(sub_contraction);
      }
      for j in i+1..vertices.len() {
        vertices[j].occupied[a].o -= assign_in[j-i-1];
      }
    }
    for j in i+1..vertices.len() {
      vertices[j].occupied[a].i -= assign_out[j-i-1];
    }
  }

  out
}

#[derive(Debug)]
pub struct Graph {
  // bond coloring
  pub adj: Vec<Vec<usize>>,
  // bond colors to degrees
  pub edge_kinds: Vec<Vec<(usize,usize)>>,
  // vert coloring
  pub nodes: Vec<usize>,
  // raw contractions
  pub contraction: Option<Contraction>,
}
impl Hash for Graph {
  fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
    let mut key = build_node_contexts(self);
    key.sort_unstable();
    for item in key {
      item.iter().for_each(|x| x.hash(state));
    }
  }
}
impl PartialEq for Graph {
  fn eq(&self, other: &Self) -> bool {
    vf2_graph_iso(self, other)
  }
}
impl Eq for Graph {}

fn build_node_contexts(g: &Graph) -> Vec<Vec<(usize,usize)>> {
  let mut key = vec![];
  for i in 0..g.nodes.len() {
    let v = g.nodes[i];
    let mut ec: Vec<(usize,usize)> = g.adj[i].iter()
      .filter(|&&x| x != 0).map(|&x| &g.edge_kinds[x]).flatten().cloned().collect();
    ec.sort_unstable();
    // HACK: add node color
    ec.push((v,0));
    key.push(ec);
  }
  key
}

fn vf2_graph_iso_core(
  g1: &Graph, g2: &Graph, key1: &Vec<Vec<(usize,usize)>>, key2: &Vec<Vec<(usize,usize)>>,
  partial_fwd: &mut Vec<Option<usize>>,
  partial_bwd: &mut Vec<Option<usize>>, next: usize,
) -> bool {
  if next == g1.nodes.len() {
    return true;
  }
  let key1_n1 = &key1[next];
  for next2 in 0..g2.nodes.len() {
    // already matched
    if partial_bwd[next2].is_some() {
      continue;
    }
    // key mismatch (incl. node color mismatch)
    let key2_n2 = &key2[next2];
    if key1_n1 != key2_n2 {
      continue;
    }
    // edge mismatch within current partial matching
    let mut edge_match = true;
    for (i,color) in g1.adj[next].iter().enumerate() {
      if i >= next {
        break;
      }
      assert!(partial_fwd[i].is_some());
      let fwd_i = partial_fwd[i].unwrap();
      if g2.edge_kinds[g2.adj[next2][fwd_i]] != g1.edge_kinds[*color] {
        edge_match = false;
        break;
      }
    }
    if !edge_match {
      continue;
    }
    partial_fwd[next] = Some(next2);
    partial_bwd[next2] = Some(next);
    let sub_match = vf2_graph_iso_core(
      g1, g2, key1, key2, partial_fwd, partial_bwd, next+1);
    if sub_match {
      return true;
    }
    partial_fwd[next] = None;
    partial_bwd[next2] = None;
  }
  return false;
}

fn vf2_graph_iso(g1: &Graph, g2: &Graph) -> bool {
  if g1.nodes.len() != g2.nodes.len() {
    return false;
  }
  // NOTE: For charged particles, care is needed.
  let key1 = build_node_contexts(g1);
  let key2 = build_node_contexts(g2);
  // for all i < next:
  // g1.nodes[i] ~ g2.nodes[partial_fwd[i]]
  // and partial_bwd[partial_fwd[i]] = i
  let mut partial_fwd = vec![None; g1.nodes.len()];
  let mut partial_bwd = vec![None; g1.nodes.len()];
  vf2_graph_iso_core(g1, g2, &key1, &key2, &mut partial_fwd, &mut partial_bwd, 0)
}

fn contraction_to_graph(
  vert_kinds: &Vec<usize>, c: &Contraction, flavors: &Vec<Flavor>
) -> Graph
{
  let n_flavor = flavors.len();
  let n_nodes = vert_kinds.len();
  let mut degrees = vec![vec![vec![Degree::new(); n_nodes]; n_nodes]; n_flavor];
  for bond in &c.bonds {
    if flavors[bond.flavor].charged {
      degrees[bond.flavor][bond.i][bond.j] = bond.degree;
      degrees[bond.flavor][bond.j][bond.i] = bond.degree.conj();
    }
    else {
      degrees[bond.flavor][bond.i][bond.j] = bond.degree;
      degrees[bond.flavor][bond.j][bond.i] = bond.degree;
    }
  }
  let mut adj = vec![vec![0; n_nodes]; n_nodes];
  let mut edge_colors = HashMap::new();
  let mut edge_kinds = vec![];
  for i in 0..n_nodes {
    for j in 0..n_nodes {
      let bond: Vec<(usize,usize)> = degrees.iter().map(
        |degree_flav| {
          let degree = degree_flav[i][j];
          (degree.i, degree.o)
        }).collect();
      if !edge_colors.contains_key(&bond) {
        let color = edge_colors.len();
        edge_kinds.push(bond.clone());
        edge_colors.insert(bond, color);
        adj[i][j] = color;
      }
      else {
        adj[i][j] = *edge_colors.get(&bond).unwrap();
      }
    }
  }
  Graph {
    adj, edge_kinds,
    nodes: vert_kinds.clone(),
    contraction: None,
  }
}

#[derive(Debug)]
pub struct SymmFactor {
  // edge symmetry factor (divide)
  pub edge_symm: BigUint,
  // vertex symmetry factor (divide)
  pub vert_symm: BigUint,
  // graph iso count (multiply)
  pub count: BigUint,
}
impl SymmFactor {
  fn new() -> Self {
    Self {
      edge_symm: BigUint::from(1u64),
      vert_symm: BigUint::from(1u64),
      count: BigUint::from(1u64),
    }
  }
  /// total symmetry factor to divide by
  pub fn total(&self) -> BigUint {
    let symm = &self.vert_symm * &self.edge_symm;
    assert!(symm.is_multiple_of(&self.count));
    symm / &self.count
  }
  /// total symmetry factor to divide by, cast to u64 or raise error
  pub fn total_u64(&self) -> Result<u64, Error> {
    let digits = self.total().to_u64_digits();
    if digits.len() == 1{
      Ok(digits[0])
    }
    else {
      Err(Error::General("numeric overflow in symmetry factor".into()))
    }
  }
}

pub fn is_connected(g: &Graph) -> bool {
  let n_nodes = g.adj.len();
  if n_nodes == 0 {
    return true;
  }
  let is_connecting: Vec<bool> = g.edge_kinds.iter()
    .map(|degrees| {
      for &(deg_in, deg_out) in degrees.iter() {
        if deg_in > 0 || deg_out > 0 {
          return true;
        }
      }
      return false;
    })
    .collect();
  let mut queue = vec![0];
  let mut seen_count = 1;
  let mut seen = vec![false; n_nodes];
  seen[0] = true;
  while queue.len() > 0 {
    let node = queue.pop().unwrap();
    for i in 0..g.adj[node].len() {
      if seen[i] {
        continue;
      }
      if !is_connecting[g.adj[node][i]] {
        continue;
      }
      seen[i] = true;
      seen_count += 1;
      queue.push(i);
    }
  }
  seen_count == n_nodes
}

type ContractionFilter = Box<dyn Fn(&Contraction, &Vec<usize>, &Theory) -> bool>;
type GraphFilter = Box<dyn Fn(&Graph) -> bool>;
pub struct Filters {
  contraction_filter: Option<ContractionFilter>,
  graph_filter: Option<GraphFilter>
}
impl Filters {
  /// No active filters
  pub fn none() -> Self {
    Self {
      contraction_filter: None,
      graph_filter: None,
    }
  }
  pub fn filter_contraction(f: ContractionFilter) -> Self {
    Self {
      contraction_filter: Some(f),
      graph_filter: None,
    }
  }
  pub fn filter_graph(f: GraphFilter) -> Self {
    Self {
      contraction_filter: None,
      graph_filter: Some(f),
    }
  }
  pub fn filter_both(cf: ContractionFilter, gf: GraphFilter) -> Self {
    Self {
      contraction_filter: Some(cf),
      graph_filter: Some(gf),
    }
  }
}

pub fn enumerate_distinct_graphs(
  vert_kinds: Vec<usize>, theory: &Theory, filters: &Filters,
) -> HashMap<Graph, SymmFactor> {
  let mut vertices = vert_kinds.iter().map(|&ind| theory.make_vertex(ind)).collect();
  let contractions = enumerate_contractions(&mut vertices, 0, 0, theory);
  let mut graphs: HashMap<Graph, SymmFactor> = HashMap::new();
  contractions.into_iter()
    .filter(|c| match &filters.contraction_filter {
      Some(f) => f(c, &vert_kinds, theory),
      None => true,
    })
    .for_each(|c| {
      let mut graph = contraction_to_graph(&vert_kinds, &c, &theory.flavors);
      if !graphs.contains_key(&graph) {
        graph.contraction = Some(c);
        graphs.insert(graph, SymmFactor::new());
      }
      else {
        graphs.get_mut(&graph).unwrap().count += BigUint::from(1u64);
      }
    });
  let mut graphs = match &filters.graph_filter {
    Some(f) => graphs.into_iter().filter(|(k,_)| f(k)).collect(),
    None => graphs,
  };

  // fill in symmetry factors
  let mut vert_symm = BigUint::from(1u64);
  let mut vert_counts = vec![0; theory.vertices.len()];
  vert_kinds.iter().for_each(|&kind| {
    vert_counts[kind] += 1;
  });
  vert_counts.iter().for_each(|&count| {
    for n in 2..=count {
      let big_n = n.to_biguint().unwrap();
      vert_symm *= big_n;
    }
  });

  graphs.iter_mut().for_each(|graph_count| {
    for bond in &graph_count.0.contraction.as_ref().unwrap().bonds {
      for n in 2..=bond.degree.i {
        let big_n = n.to_biguint().unwrap();
        graph_count.1.edge_symm *= big_n;
      }
      for n in 2..=bond.degree.o {
        let big_n = n.to_biguint().unwrap();
        graph_count.1.edge_symm *= big_n;
      }
    }
    graph_count.1.vert_symm = vert_symm.clone();
  });
  graphs
}

#[cfg(test)]
mod tests {
  use super::*;
  use std::collections::hash_map::DefaultHasher;

  #[test]
  fn partition_into_simple() {
    let n = 10;
    let sizes: Vec<usize> = vec![3, 4, 5, 6];
    let parts = partition_into(n, &sizes);
    for part in parts.iter() {
      assert_eq!(part.iter().sum::<usize>(), n);
      for i in 0..part.len() {
        assert!(part[i] <= sizes[i]);
      }
    }
    assert_eq!(parts.len(), 96);
  }

  #[test]
  fn partition_into_zeros() {
    let n = 0;
    let sizes: Vec<usize> = vec![0, 0, 0, 0];
    let parts = partition_into(n, &sizes);
    assert_eq!(parts, vec![vec![0, 0, 0, 0]]);
  }

  #[test]
  fn real_scalar() {
    let real_scalar = Flavor { stats: FermiStats::Boson, charged: false };
    let flavors = vec![real_scalar];
    let cubic = VertexKind { degrees: vec![Degree { i: 3, o: 0 }] };
    let vertices = vec![cubic];
    let theory = Theory {
      flavors, vertices
    };

    let mut one = vec![theory.make_vertex(0)];
    let res1 = enumerate_contractions(&mut one, 0, 0, &theory);
    assert_eq!(res1.len(), 0);

    let mut two = vec![
      theory.make_vertex(0),
      theory.make_vertex(0),
    ];
    let res2 = enumerate_contractions(&mut two, 0, 0, &theory);
    assert_eq!(res2.len(), 1);
    let contraction = &res2[0];
    assert_eq!(contraction.bonds.len(), 1);
    let bond = &contraction.bonds[0];
    assert_eq!(bond, &Bond {
      flavor: 0, degree: Degree { i: 3, o: 0 },
      i: 0, j: 1
    });

    let mut three = vec![
      theory.make_vertex(0),
      theory.make_vertex(0),
      theory.make_vertex(0),
    ];
    let res3 = enumerate_contractions(&mut three, 0, 0, &theory);
    assert_eq!(res3.len(), 0);

    let mut four = vec![
      theory.make_vertex(0),
      theory.make_vertex(0),
      theory.make_vertex(0),
      theory.make_vertex(0),
    ];
    let res4 = enumerate_contractions(&mut four, 0, 0, &theory);
    assert_eq!(res4.len(), 10);

    let graphs4 = enumerate_distinct_graphs([0,0,0,0].to_vec(), &theory, &Filters::none());
    assert_eq!(graphs4.len(), 3);
    // TODO: more stringent checks
  }

  #[test]
  fn vf2_six_node_case1() {
    let g1 = Graph {
      nodes: vec![0; 6],
      adj: vec![
        vec![0,2,0,0,0,0],
        vec![2,0,0,0,0,0],
        vec![0,0,0,1,1,1],
        vec![0,0,1,0,1,1],
        vec![0,0,1,1,0,1],
        vec![0,0,1,1,1,0],
      ],
      edge_kinds: vec![
        vec![(0,0)],
        vec![(1,0)],
        vec![(3,0)],
      ],
      contraction: None,
    };
    let g2 = Graph {
      nodes: vec![0; 6],
      adj: vec![
        vec![0,0,0,2,0,0],
        vec![0,0,1,0,1,1],
        vec![0,1,0,0,1,1],
        vec![2,0,0,0,0,0],
        vec![0,1,1,0,0,1],
        vec![0,1,1,0,1,0],
      ],
      edge_kinds: vec![
        vec![(0,0)],
        vec![(1,0)],
        vec![(3,0)],
      ],
      contraction: None,
    };
    let g3 = Graph {
      nodes: vec![0; 6],
      adj: vec![
        vec![0,0,2,0,2,2],
        vec![0,0,0,1,0,0],
        vec![2,0,0,0,2,2],
        vec![0,1,0,0,0,0],
        vec![2,0,2,0,0,2],
        vec![2,0,2,0,2,0],
      ],
      edge_kinds: vec![
        vec![(0,0)],
        vec![(3,0)],
        vec![(1,0)],
      ],
      contraction: None,
    };

    assert_eq!(g1, g1);
    assert_eq!(g2, g2);
    assert_eq!(g3, g3);
    assert_eq!(g1, g2);
    assert_eq!(g1, g3);
    assert_eq!(g2, g3);

    let mut hasher = DefaultHasher::new();
    let hash1 = g1.hash(&mut hasher);
    let mut hasher = DefaultHasher::new();
    let hash2 = g2.hash(&mut hasher);
    assert_eq!(hash1, hash2);
  }

  #[test]
  fn two_to_two_scalar_tree_level() {
    let real_scalar = Flavor { stats: FermiStats::Boson, charged: false };
    let flavors = vec![real_scalar];
    let cubic = VertexKind { degrees: vec![Degree { i: 3, o: 0 }] };
    let external = VertexKind { degrees: vec![Degree { i: 1, o: 0 }] };
    let vertices = vec![cubic, external];
    let theory = Theory {
      flavors, vertices
    };

    let graphs = enumerate_distinct_graphs(vec![0,0,1,1,1,1], &theory, &Filters::none());
    // two disconnected and one connected topology
    assert_eq!(graphs.len(), 3);
  }

  #[test]
  fn two_to_two_qed_tree_level() {
    let scalar = Flavor { stats: FermiStats::Boson, charged: true };
    let photon = Flavor { stats: FermiStats::Boson, charged: false };
    let flavors = vec![scalar, photon];
    let interaction = VertexKind { degrees: vec![
      Degree { i: 1, o: 1 }, // matter
      Degree { i: 1, o: 0 }, // photon
    ]};
    let ext_in = VertexKind { degrees: vec![
      Degree { i: 0, o: 1 }, Degree::new()
    ]};
    let ext_out = VertexKind { degrees: vec![
      Degree { i: 1, o: 0 }, Degree::new()
    ]};
    let vertices = vec![interaction, ext_in, ext_out];
    let theory = Theory {
      flavors, vertices
    };

    let graphs = enumerate_distinct_graphs(vec![0,0,1,1,2,2], &theory, &Filters::none());
    // connected + singly disconnected + doubly disconnected
    assert_eq!(graphs.len(), 3);
  }

  #[test]
  fn two_to_two_qed_tree_level_conn() {
    let scalar = Flavor { stats: FermiStats::Boson, charged: true };
    let photon = Flavor { stats: FermiStats::Boson, charged: false };
    let flavors = vec![scalar, photon];
    let interaction = VertexKind { degrees: vec![
      Degree { i: 1, o: 1 }, // matter
      Degree { i: 1, o: 0 }, // photon
    ]};
    let ext_in = VertexKind { degrees: vec![
      Degree { i: 0, o: 1 }, Degree::new()
    ]};
    let ext_out = VertexKind { degrees: vec![
      Degree { i: 1, o: 0 }, Degree::new()
    ]};
    let vertices = vec![interaction, ext_in, ext_out];
    let theory = Theory {
      flavors, vertices
    };

    let graphs = enumerate_distinct_graphs(
      vec![0,0,1,1,2,2], &theory,
      &Filters::filter_graph(Box::new(is_connected)));
    // only one connected diagram up to iso
    assert_eq!(graphs.len(), 1);
  }
}
