User guide
==========

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

        .. button-ref:: boolean_operations
            :expand:
            :color: primary
            :click-parent:

            To Boolean ops

    .. grid-item-card::
        :img-top: ../source/_static/index-images/step.svg


        ^^^^^^^^^^^^^^

        Dealing with STEP files

        +++

        .. button-ref:: dealing_with_step
            :expand:
            :color: primary
            :click-parent:

            STEP files

    .. grid-item-card::
        :img-top: ../source/_static/index-images/stl_file.svg


        ^^^^^^^^^^^^^^

        Dealing with STL files

        +++

        .. button-ref:: dealing_with_stl
            :expand:
            :color: primary
            :click-parent:

            STL files

    .. grid-item-card::
        :img-top: ../source/_static/index-images/tools-folder-svgrepo-com.svg


        ^^^^^^^^^^^^^^

        Very useful Operations

        +++

        .. button-ref:: useful_operations
            :expand:
            :color: primary
            :click-parent:

            Useful

The volmdlr library is structured using a modular approach as follows:


* :ref:`core_compiled`: provides fundamental objects and algorithms for B-Rep modeling, implemented in C++ for performance reasons;
* :ref:`core`: provides fundamental data structures and algorithms for B-Rep modeling;
* :ref:`geometry`: provides functions and tools for calculations and geometric operations with volmdlr geometric objects;
* :ref:`curves`: Provides fundamental curves Objects, like Infinite lines, Circle and ellipses.
* :ref:`edges`: a shape corresponding to a curve and bounded by a start and an end point;
* :ref:`wires`: a sequence of edges connected by their vertices;
* :ref:`surfaces`: Provides surfaces objects, which are defined by its mathematical equations
* :ref:`faces`: Provides faces objects, whire are defined by a surface and a boundary representation, that is called Suface2D in volmdlr;
* :ref:`shells`: Provides shells objects, whiche are collections of faces
* :ref:`stl`: provides support for importing and exporting STL files;
* :ref:`step`: provides support for importing and exporting STEP files;
* :ref:`mesh`: provides tools for creating and manipulating triangular meshes;
* :ref:`primitives3d`: provides tools for creating 3D primitives like extrusion, revolution and sweep.
* :ref:`display`: provides tools for rendering and visualizing 3D models;
* :ref:`cloud`: provides algorithms to work with a cloud of points;


:ref:`modules`
--------------


Q&A
---

No questions yet.
