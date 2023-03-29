Getting started
===============

Install
-------

.. code::

  pip(3) install Volmdlr

Important
---------

Before using Volmdlr, be sure to have a compiler C/C++ (not necessary on Linux).
N.B : With Windows you have to download one and allows it to read Python's code.

Introduction
---------
The volmdlr library is an open-source Python library primarily developed by Dessia Technologies,
aimed at providing 3D modeling capabilities based on Boundary Representation (B-Rep) purely in Python.
The library is designed to be easy to use, efficient, and customizable.

Boundary Representation is a technique used in computer-aided design (CAD) and computer-aided
manufacturing (CAM) systems to represent the geometry of a solid object as a collection of surfaces
and curves. B-Rep is widely used in CAD and CAM due to its ability to accurately represent complex
shapes and to provide a rich set of operations for manipulating and analyzing these shapes.

The volmdlr library is structured using a modular approach as follows:

* **cloud**: provides algorithms to work with a cloud of points.
* **core_compiled**: provides fundamental objects and algorithms for B-Rep modeling, implemented in C++ for performance reasons
* **core**: provides fundamental data structures and algorithms for B-Rep modeling;
* **display**: provides tools for rendering and visualizing 3D models;
* **geometry**: provides functions and tools for calculations and geometric operations with volmdlr geometric objects;
* **edges**: a shape corresponding to a curve and bounded by a start and an end point;
* **wires**:a sequence of edges connected by their vertices;
* **faces**: This module regroup surfaces, faces and shell objects. It's supposed to be devided in three different modules in future releases. A surface is defined by its mathematical equations. A face is defined by a surface and a boundary representation, that is called Suface2D in volmdlr. And a Shell is a collection of faces;
* **stl**: provides support for importing and exporting STL files;
* **step**: provides support for importing and exporting STEP files;
* **mesh**: provides tools for creating and manipulating triangular meshes;
* **primitives3d**: provides tools for creating 3D primitives like extrusion, revolution and sweep;


CAD 
---

The following picture is here to show you about the power of Volmdlr.

.. image:: images/Castest.jpg
  :width: 600

https://github.com/Dessia-tech/volmdlr/blob/distancewire/scripts/Tutorial/CasTest.py

This engine is made from ExtrudedProfile and RoundedLineSegment2D. With CylindricalFace3D and ToroidalFace3D too.
You can see Sweep as presented before.
If you want to load it, run the code above, and if you are curious, check on Primitives3D how to create Faces.
