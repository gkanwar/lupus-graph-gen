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
  let mut four = vec![
    theory.make_vertex(0),
    theory.make_vertex(0),
    theory.make_vertex(0),
    theory.make_vertex(0),
  ];
  let res4 = enumerate_contractions(&mut four, 0, 0, &theory);
  let f = File::create("test.pdf").map_err(Error::IOError)?;
  let mut pdf = PdfWriter::new(f)?;
  for contraction in &res4 {
    draw_contraction(contraction, &four, &theory, &mut pdf)?;
    println!("{:?}", contraction);
  }
  pdf.write_all()?;
  Ok(())
}
