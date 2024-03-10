Lupus graph generator
=====================
A stab at generating Feynman graphs in a general and (hopfully!) readable way
using a very direct Rust implementation. This is a pedagogical exercise and not
intended for usage in real work.

The basic algorithm is similar to prior work:
 * All graphs over an ordered list of nodes is generated
 * Isomorphic graphs are collected together
 * A total symmetry factor can be obtained from the per-bond symmetries
   and the count of isomorphic graphs


Features
========
 - [x] Automatic PDF outputs
 - [x] Multiple flavors (needs testing)
 - [x] Charged and uncharged (needs testing)
 - [ ] Fermions
 - [x] Vertex and edge symmetry factors
 - [ ] Tadpoles

