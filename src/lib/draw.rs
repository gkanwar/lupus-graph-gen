/// Module to draw Feynman graphs as PDF content streams.

use crate::core::{Graph, Theory};
use crate::error::Error;
use std::io::{self, Seek, Write};
use std::result::Result;
use std::f64::consts::PI;

// some shorthands
type Res = Result<(), Error>;
const OK: Res = Ok(());

struct PdfObject {
  ind: usize,
  content: PdfValue,
}
enum PdfValue {
  Name(String),
  Indirect(usize),
  Boolean(bool),
  Integer(i64),
  Real(f64),
  Array(Vec<PdfValue>),
  Dict(Vec<(String, PdfValue)>),
  ContentStream(String),
}
impl From<i64> for PdfValue {
  fn from(value: i64) -> Self {
    Self::Integer(value)
  }
}
impl<T: Into<PdfValue>> From<Vec<T>> for PdfValue
where PdfValue: From<T>
{
  fn from(values: Vec<T>) -> Self {
    PdfValue::Array(values.into_iter().map(|v| PdfValue::from(v)).collect())
  }
}
impl ToString for PdfValue {
  fn to_string(&self) -> String {
    match self {
      Self::Name(name) => format!("/{}", name),
      Self::Indirect(obj) => format!("{} 0 R", obj+1),
      Self::Boolean(value) => match value {
        true => "true".into(),
        false => "false".into(),
      },
      Self::Integer(value) => format!("{}", value),
      Self::Real(value) => format!("{}", value),
      Self::Array(values) => {
        let mut out = String::new();
        out += "[ ";
        for value in values {
          out += &value.to_string();
          out += " ";
        }
        out += "]";
        out
      }
      Self::Dict(items) => {
        let mut out = String::new();
        out += "<<\r\n";
        for (k,v) in items {
          assert_eq!(k.chars().nth(0).unwrap(), '/');
          out += &format!("{} {}\r\n", k, v.to_string());
        }
        out += ">>";
        out
      }
      Self::ContentStream(s) => {
        let mut out = String::new();
        out += &PdfValue::Dict(vec![
          ("/Length".into(), (s.len() as i64).into())
        ]).to_string();
        out += "\r\nstream\r\n";
        out += &s;
        out += "\r\nendstream";
        out
      }
    }
  }
}

fn write_object(obj: &PdfObject, f: &mut impl Write) -> Res {
  f.write(format!("{} 0 obj\r\n", obj.ind+1).as_bytes()).map_err(Error::IOError)?;
  f.write(obj.content.to_string().as_bytes()).map_err(Error::IOError)?;
  f.write("\r\nendobj\r\n".as_bytes()).map_err(Error::IOError)?;
  return OK;
}


pub struct PdfWriter<F: Seek + Write> {
  objects: Vec<PdfObject>,
  root_obj: usize,
  info_obj: usize,
  pages_obj: usize,
  font_objs: Vec<usize>,
  obj_offsets: Option<Vec<u64>>,
  xref_offset: Option<u64>,
  f: F,
}

impl<F: Seek + Write> PdfWriter<F> {
  pub fn new(mut f: F) -> Result<PdfWriter<F>, Error> {
    f.seek(io::SeekFrom::Start(0)).map_err(Error::IOError)?;
    let mut writer = PdfWriter::<F> {
      objects: vec![],
      root_obj: 0,
      info_obj: 0,
      pages_obj: 0,
      font_objs: vec![],
      obj_offsets: None,
      xref_offset: None,
      f
    };
    writer.pages_obj = writer.create_object(PdfValue::Dict(
      vec![
        ("/Type".into(), PdfValue::Name("Pages".into())),
        ("/Kids".into(), PdfValue::Array(vec![])),
        ("/Count".into(), 0.into()),
      ]
    ));
    writer.root_obj = writer.create_object(PdfValue::Dict(
      vec![
        ("/Type".into(), PdfValue::Name("Catalog".into())),
        ("/Pages".into(), PdfValue::Indirect(writer.pages_obj)),
      ]
    ));
    writer.info_obj = writer.create_object(PdfValue::Dict(
      vec![]
    ));
    Ok(writer)
  }

  fn create_object(&mut self, content: PdfValue) -> usize {
    let ind = self.objects.len();
    self.objects.push(PdfObject {
      ind, content
    });
    ind
  }

  fn create_page(&mut self, content: PdfValue, bounds: [i64; 4]) -> Res {
    // assert!(match resources { PdfValue::Dict(_) => true, _ => false });
    let mut resources = vec![];
    for i in 0..self.font_objs.len() {
      resources.push(
        (format!("/F{}", i+1), PdfValue::Indirect(self.font_objs[i])));
    }
    let page = PdfValue::Dict(
      vec![
        ("/Type".into(), PdfValue::Name("Page".into())),
        ("/MediaBox".into(), bounds.into_iter().collect::<Vec<i64>>().into()),
        ("/Resources".into(), PdfValue::Dict(resources)),
        ("/Parent".into(), PdfValue::Indirect(self.pages_obj)),
        ("/Contents".into(), content),
      ]
    );
    let ind = self.create_object(page);
    let pages = &mut self.objects[self.pages_obj];
    match &mut pages.content {
      PdfValue::Dict(items) => {
        match &mut items[1] {
          (key, PdfValue::Array(pages)) => {
            assert_eq!(key, "/Kids");
            pages.push(PdfValue::Indirect(ind))
          },
          _ => { unreachable!(); }
        }
        match &mut items[2] {
          (key, PdfValue::Integer(count)) => {
            assert_eq!(key, "/Count");
            *count += 1;
          }
          _ => { unreachable!(); }
        }
      }
      _ => { unreachable!(); }
    }
    return OK;
  }

  pub fn add_font(&mut self, name: String) -> Res {
    let font_dict = PdfValue::Dict(vec![
      ("/Type".into(), PdfValue::Name("Font".into())),
      ("/Subtype".into(), PdfValue::Name("Type1".into())),
      ("/BaseFont".into(), PdfValue::Name(name)),
    ]);
    let obj = self.create_object(font_dict);
    self.font_objs.push(obj);
    return OK;
  }

  pub fn write_all(&mut self) -> Res {
    self.write_header()?;
    self.write_objects()?;
    self.write_xref_table()?;
    self.write_trailer()?;
    return OK;
  }

  fn write_header(&mut self) -> Res {
    // PDF version
    self.f.write("%PDF-1.4\r\n".as_bytes()).map_err(Error::IOError)?;
    // PDF spec requires a comment line with arbitary binary to ensure file inspectors
    // correctly guess that this file contains binary content
    self.f.write(&[b'%', 0xff, 0xff, 0xff, 0xff, b'\r', b'\n']).map_err(Error::IOError)?;
    return OK;
  }

  fn write_objects(&mut self) -> Res {
    let mut offsets = vec![];
    for i in 0..self.objects.len() {
      offsets.push(self.f.stream_position().map_err(Error::IOError)?);
      write_object(&self.objects[i], &mut self.f)?;
    }
    self.obj_offsets = Some(offsets);
    return OK;
  }

  fn write_xref_table(&mut self) -> Res {
    assert!(self.xref_offset.is_none());
    self.xref_offset = Some(self.f.stream_position().map_err(Error::IOError)?);
    self.f.write("xref\r\n".as_bytes()).map_err(Error::IOError)?;
    let n_obj = self.obj_offsets.as_ref().unwrap().len()+1;
    self.f.write(format!("0 {}\r\n", n_obj).as_bytes()).map_err(Error::IOError)?;
    self.f.write("0000000000 65535 f\r\n".as_bytes()).map_err(Error::IOError)?;
    for off in self.obj_offsets.as_ref().unwrap() {
      self.f.write(format!("{:010} 00000 n\r\n", off).as_bytes()).map_err(Error::IOError)?;
    }
    return OK;
  }

  fn write_trailer(&mut self) -> Res {
    let n_obj = self.objects.len()+1;
    let xref_offset = self.xref_offset.unwrap();
    // TODO: use PdfValues
    self.f.write("trailer\r\n<<\r\n".as_bytes()).map_err(Error::IOError)?;
    self.f.write(format!("/Size {}\r\n", n_obj).as_bytes()).map_err(Error::IOError)?;
    self.f.write(format!("/Root {} 0 R\r\n", self.root_obj+1).as_bytes()).map_err(Error::IOError)?;
    self.f.write(format!("/Info {} 0 R\r\n", self.info_obj+1).as_bytes()).map_err(Error::IOError)?;
    // self.f.write("/ID [ <deadbeef> ] [ <deadbeef> ]\r\n".as_bytes()).map_err(Error::IOError)?;
    self.f.write(">>\r\n".as_bytes()).map_err(Error::IOError)?;
    self.f.write("startxref\r\n".as_bytes()).map_err(Error::IOError)?;
    self.f.write(format!("{}\r\n", xref_offset).as_bytes()).map_err(Error::IOError)?;
    self.f.write("%%EOF\r\n".as_bytes()).map_err(Error::IOError)?;
    return OK;
  }
}

type Coord = (f64,f64);

/// Convert drawing coords to coords printer points
fn to_points(coord: Coord) -> Coord {
  let x = coord.0;
  let y = coord.1;
  (x*72.0, y*72.0)
}

fn push_coord(coord: Coord, stream: &mut String) {
  let pt = to_points(coord);
  stream.push_str(&format!("{} {} ", pt.0, pt.1));
}

fn draw_node(coord: Coord, stream: &mut String) {
  const DX: f64 = 0.1;
  push_coord((coord.0-DX, coord.1-DX), stream);
  stream.push_str("m ");
  push_coord((coord.0-DX, coord.1+DX), stream);
  stream.push_str("l ");
  push_coord((coord.0+DX, coord.1+DX), stream);
  stream.push_str("l ");
  push_coord((coord.0+DX, coord.1-DX), stream);
  stream.push_str("l h ");

  stream.push_str("f ");
}

fn draw_line(x: Coord, y: Coord, stream: &mut String) {
  push_coord(x, stream);
  stream.push_str("m ");
  push_coord(y, stream);
  stream.push_str("l S ");
}

fn draw_bond(x: Coord, y: Coord, degree: usize, stream: &mut String) {
  if degree == 1 {
    draw_line(x, y, stream);
    return;
  }
  let mut u: Coord = (y.0 - x.0, y.1 - x.1);
  let norm = (u.0*u.0 + u.1*u.1).sqrt();
  u.0 /= norm;
  u.1 /= norm;
  let v: Coord = (u.1, -u.0);
  for i in 0..degree {
    let mut r = (i as f64) / ((degree-1) as f64); // in [0,1]
    r -= 0.5;
    r *= 0.1;
    let xp = (x.0 + r*v.0, x.1 + r*v.1);
    let yp = (y.0 + r*v.0, y.1 + r*v.1);
    draw_line(xp, yp, stream);
  }
}

fn config_drawing(stream: &mut String) {
  stream.push_str("3 w 1 J 1 j ");
}


pub fn draw_contraction<F: Seek + Write>(
  graph: &Graph, count: usize,
  theory: &Theory, pdf: &mut PdfWriter<F>
) -> Res {
  let c = &graph.contraction.as_ref().unwrap();
  // TODO: externals
  // FORNOW: assume vacuum bubble
  // FORNOW: construct a circle geometry just so we see everything
  let mut stream = String::new();
  config_drawing(&mut stream);

  // draw nodes
  let mut nodes = vec![];
  for i in 0..graph.nodes.len() {
    // TODO: do something with vertex type
    let v = &graph.nodes[i];
    let theta = 2.0*PI*(i as f64) / (graph.nodes.len() as f64);
    nodes.push((theta.cos()+1.5, theta.sin()+1.5));
    draw_node(nodes[nodes.len()-1], &mut stream);
  }

  // draw bonds
  for bond in &c.bonds {
    if bond.degree.i == 0 && bond.degree.o == 0 {
      continue;
    }
    let x = nodes[bond.i];
    let y = nodes[bond.j];
    // TODO: flavors
    // TODO: better degree handling
    draw_bond(x, y, bond.degree.i+bond.degree.o, &mut stream);
  }

  // Note: Assumes at least one font has been pushed
  stream.push_str("BT /F1 12 Tf ");
  push_coord((0.1, 0.1), &mut stream);
  stream.push_str("Td ");
  stream.push_str(&format!("(S~{}) Tj ", count));
  stream.push_str("ET ");
  
  let content = PdfValue::Indirect(pdf.create_object(PdfValue::ContentStream(stream)));
  pdf.create_page(content, [0, 0, 72*3, 72*3])?;
  
  return OK;
}
