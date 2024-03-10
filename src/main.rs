use lupus::prelude::*;
use std::fs::File;

fn main() -> Result<(), Error> {
  let real_scalar = Flavor { stats: FermiStats::Boson, charged: false };
  let flavors = vec![real_scalar];
  let cubic = VertexKind { degrees: vec![Degree { i: 3, o: 0 }] };
  let vertices = vec![cubic];
  let theory = Theory {
    flavors, vertices
  };
  {
    let verts = vec![0,0,0,0];
    let graphs = enumerate_distinct_graphs(verts, &theory);
    let f = File::create("scalar_4_phi4.pdf").map_err(Error::IOError)?;
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
  }
  {
    let verts = vec![0,0,0,0,0,0];
    let graphs = enumerate_distinct_graphs(verts, &theory);
    let f = File::create("scalar_6_phi4.pdf").map_err(Error::IOError)?;
    let mut pdf = PdfWriter::new(f)?;
    for (graph,symm) in &graphs {
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
  }
  {
    let verts = vec![0,0,0,0,0,0,0,0];
    let graphs = enumerate_distinct_graphs(verts, &theory);
    let f = File::create("scalar_8_phi4.pdf").map_err(Error::IOError)?;
    let mut pdf = PdfWriter::new(f)?;
    for (graph,symm) in &graphs {
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
  }
  Ok(())
}
