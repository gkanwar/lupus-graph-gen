use lupus::prelude::*;
use std::fs::File;

fn run_case(fname: String, verts: Vec<usize>, theory: &Theory) -> Result<(), Error> {
  let graphs = enumerate_distinct_graphs(verts, &theory);
  let f = File::create(fname).map_err(Error::IOError)?;
  let mut pdf = PdfWriter::new(f)?;
  for (graph,symm) in &graphs {
    println!("{:?}", symm);
    let factor = symm.total();
    let factor = if factor.to_u64_digits().len() == 1 {
      factor.to_u64_digits()[0] as usize
    }
    else {
      return Err(Error::General("Count too large".into()));
    };
    draw_contraction(graph, factor, &theory, &mut pdf)?;
    println!("{:?} (S=1/{})", graph, factor);
  }
  pdf.write_all()?;
  return Ok(())
}

fn main() -> Result<(), Error> {
  {
    let real_scalar = Flavor { stats: FermiStats::Boson, charged: false };
    let flavors = vec![real_scalar];
    let cubic = VertexKind { degrees: vec![Degree { i: 3, o: 0 }] };
    let external = VertexKind { degrees: vec![Degree { i: 1, o: 0 }] };
    let vertices = vec![cubic, external];
    let scalar_phi3 = Theory {
      flavors, vertices
    };
    run_case("scalar_2to2_phi3.pdf".into(), vec![0,0,1,1,1,1], &scalar_phi3)?;
    run_case("scalar_4_phi3.pdf".into(), vec![0,0,0,0], &scalar_phi3)?;
    run_case("scalar_6_phi3.pdf".into(), vec![0,0,0,0,0,0], &scalar_phi3)?;
    run_case("scalar_8_phi3.pdf".into(), vec![0,0,0,0,0,0,0,0], &scalar_phi3)?;
  }
  {
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
    let qed = Theory {
      flavors, vertices
    };
    run_case("qed_2to2.pdf".into(), vec![0,0,1,1,2,2], &qed)?;
  }
  Ok(())
}
