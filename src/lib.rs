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
#[derive(Debug,Clone,PartialEq)]
pub struct Degree {
  pub i: usize,
  pub o: usize,
}
impl Degree {
  pub fn new() -> Self {
    Self { i: 0, o: 0 }
  }
}
#[derive(Debug)]
struct VertexKind {
  degrees: Vec<Degree>,
}
#[derive(Debug)]
pub struct VertexState {
  occupied: Vec<Degree>,
  kind: usize,
}

#[derive(Debug)]
pub struct Theory {
  flavors: Vec<Flavor>,
  vertices: Vec<VertexKind>,
  // TODO:
  // externals: Vec<VertexKind>,
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
  println!("args {:?} {} {}", vertices, i, a);
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
  let open_x = kind.degrees[a].i - state.occupied[a].i;
  let open_y = kind.degrees[a].o - state.occupied[a].o;
  let mut other_open_x = vec![];
  let mut other_open_y = vec![];
  for vp in vertices[i+1..].iter() {
    // for charged wire in <- out and out -> in
    if is_charged {
      other_open_x.push(kind.degrees[a].o - vp.occupied[a].o);
      other_open_y.push(kind.degrees[a].i - vp.occupied[a].i);
    }
    // for uncharged wire in <-> in and ignore out
    // (we build it anyway for uniformity)
    else {
      other_open_x.push(kind.degrees[a].i - vp.occupied[a].i);
      assert_eq!(vp.occupied[a].o, 0);
      assert_eq!(kind.degrees[a].o, 0);
      other_open_y.push(0);
    }
  }
  println!("Assigning for {} {} -> {:?} {:?}", open_x, open_y, other_open_x, other_open_y);

  let assignments_x = partition_into(open_x, &other_open_x);
  let assignments_y = partition_into(open_y, &other_open_y);
  let mut out = vec![];
  for assign_x in assignments_x.iter() {
    assert_eq!(assign_x.len(), vertices.len() - i - 1);
    for j in i+1..vertices.len() {
      vertices[j].occupied[a].i += assign_x[j-i-1];
    }
    for assign_y in assignments_y.iter() {
      assert_eq!(assign_y.len(), vertices.len() - i - 1);
      for j in i+1..vertices.len() {
        vertices[j].occupied[a].o += assign_y[j-i-1];
      }
      let mut contraction_i = vec![];
      for j in i+1..vertices.len() {
        contraction_i.push(Bond {
          flavor: a,
          degree: Degree { i: assign_x[j-i-1], o: assign_y[j-i-1] },
          i, j
        });
      }
      println!("Recursing... ");
      let sub = enumerate_contractions(vertices, i, a+1, theory);
      println!("... result {:?}", sub);
      for mut sub_contraction in sub {
        sub_contraction.bonds.extend(contraction_i.clone());
        out.push(sub_contraction);
      }
      for j in i+1..vertices.len() {
        vertices[j].occupied[a].o -= assign_y[j-i-1];
      }
    }
    for j in i+1..vertices.len() {
      vertices[j].occupied[a].i -= assign_x[j-i-1];
    }
  }
  
  out
}

#[cfg(test)]
mod tests {
  use super::*;
  #[test]
  fn partition_into_simple() {
    let n = 10;
    let sizes: Vec<usize> = vec![3, 4, 5, 6];
    let parts = partition_into(n, &sizes);
    for part in parts.iter() {
      println!("{:?}", part);
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
    // TODO: is this correct?
    assert_eq!(res4.len(), 10);
  }
}
