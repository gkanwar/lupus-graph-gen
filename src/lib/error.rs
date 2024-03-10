use std::io;

#[derive(Debug)]
pub enum Error {
  IOError(io::Error)
}
