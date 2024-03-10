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
  let four = vec![0,0,0,0];
  let graphs4 = enumerate_distinct_graphs(four, &theory);
  let f = File::create("test.pdf").map_err(Error::IOError)?;
  let mut pdf = PdfWriter::new(f)?;
  for (graph,count) in &graphs4 {
    draw_contraction(graph, *count, &theory, &mut pdf)?;
    println!("{:?} (S=1/{})", graph, count);
  }
  pdf.write_all()?;
  Ok(())
}
