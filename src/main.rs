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
    for (graph,count) in &graphs {
      draw_contraction(graph, *count, &theory, &mut pdf)?;
      println!("{:?} (S=1/{})", graph, count);
    }
    pdf.write_all()?;
  }
  {
    let verts = vec![0,0,0,0,0,0];
    let graphs = enumerate_distinct_graphs(verts, &theory);
    let f = File::create("scalar_6_phi4.pdf").map_err(Error::IOError)?;
    let mut pdf = PdfWriter::new(f)?;
    for (graph,count) in &graphs {
      draw_contraction(graph, *count, &theory, &mut pdf)?;
      println!("{:?} (S=1/{})", graph, count);
    }
    pdf.write_all()?;
  }
  {
    let verts = vec![0,0,0,0,0,0,0,0];
    let graphs = enumerate_distinct_graphs(verts, &theory);
    let f = File::create("scalar_8_phi4.pdf").map_err(Error::IOError)?;
    let mut pdf = PdfWriter::new(f)?;
    for (graph,count) in &graphs {
      draw_contraction(graph, *count, &theory, &mut pdf)?;
      println!("{:?} (S=1/{})", graph, count);
    }
    pdf.write_all()?;
  }
  Ok(())
}
