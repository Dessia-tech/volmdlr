User guide
==========

General description
-------------------
The volmdlr library is an open-source Python library primarily developed by Dessia Technologies,
aimed at providing 3D modeling capabilities based on Boundary Representation (B-Rep) purely in Python.
The library is designed to be easy to use, efficient, and customizable.

Boundary Representation is a technique used in computer-aided design (CAD) and computer-aided
manufacturing (CAM) systems to represent the geometry of a solid object as a collection of surfaces
and curves. B-Rep is widely used in CAD and CAM due to its ability to accurately represent complex
shapes and to provide a rich set of operations for manipulating and analyzing these shapes.

The volmdlr library is structured using a modular approach as follows:


* **core_compiled**: provides fundamental objects and algorithms for B-Rep modeling, implemented in C++ for performance reasons;
* **core**: provides fundamental data structures and algorithms for B-Rep modeling;
* **geometry**: provides functions and tools for calculations and geometric operations with volmdlr geometric objects;
* **curves**: Provides fundamental curves Objects, like Infinite lines, Circle and ellipses.
* **edges**: a shape corresponding to a curve and bounded by a start and an end point;
* **wires**: a sequence of edges connected by their vertices;
* **surfaces**: Provides surfaces objects, which are defined by its mathematical equations
* **faces**:. Provides faces objects, whire are defined by a surface and a boundary representation, that is called Suface2D in volmdlr;
* **shelss**: Provides shells objects, whiche are collections of faces
* **stl**: provides support for importing and exporting STL files;
* **step**: provides support for importing and exporting STEP files;
* **mesh**: provides tools for creating and manipulating triangular meshes;
* **primitives3d**: provides tools for creating 3D primitives like extrusion, revolution and sweep.
* **display**: provides tools for rendering and visualizing 3D models;
* **cloud**: provides algorithms to work with a cloud of points;

Fundamentals and usage
----------------------
.. toctree::
   :caption: Volmdlr
   :maxdepth: 2

   volmdlr

Q&A
---

No questions yet.
