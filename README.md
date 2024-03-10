Lupus graph generator
=====================
A stab at generating Feynman graphs is a general and readable(!) way using
a simplistic Rust implementation. Very much a pedagogical exercise and not
intended for usage in real work.

The basic algorithm is similar to prior work:
 * All graphs over an ordered list of nodes is generated
 * Isomorphic graphs are collected together
 * A total symmetry factor can be obtained from the per-bond symmetries
   and the count of isomorphic graphs
