User guide
==========

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

.. grid:: 3

    .. grid-item-card::
        :img-top: ../source/_static/index-images/geometric-svgrepo-com.svg

        ^^^^^^^^^^^^^^

        How to

        +++

        .. button-ref:: how-tos
            :expand:
            :color: primary
            :click-parent:

            How to

    .. grid-item-card::
        :img-top: ../source/_static/index-images/distance.svg


        ^^^^^^^^^^^^^^

        Distances operations

        +++

        .. button-ref:: distances
            :expand:
            :color: primary
            :click-parent:

            To distances

    .. grid-item-card::
        :img-top: ../source/_static/index-images/intersection-svgrepo-com.svg


        ^^^^^^^^^^^^^^

        Intersections

        +++

        .. button-ref:: intersections
            :expand:
            :color: primary
            :click-parent:

            To intersections

    .. grid-item-card::
        :img-top: ../source/_static/index-images/union-svgrepo-com.svg


        ^^^^^^^^^^^^^^

        Boolean Operations

        +++

        .. button-ref:: intersections
            :expand:
            :color: primary
            :click-parent:

            To Boolean ops

Fundamentals and usage
----------------------

.. toctree::
   :caption: Volmdlr Basics
   :maxdepth: 2

   volmdlr.utils
   modules

Q&A
---

No questions yet.
