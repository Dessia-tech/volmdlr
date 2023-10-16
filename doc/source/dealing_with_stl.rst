===
STL
===

An STL file, short for "stereolithography" or "standard tessellation language," is a widely used file format in 3D printing and computer-aided design (CAD).
It is used to represent the surface geometry of a 3D object as a collection of interconnected triangles.
This format is essential for 3D printers and other software to understand and reproduce the physical object accurately.

An STL file represents a 3D object as a mesh composed of triangular facets. These triangles collectively define the surface geometry of the object.
Each triangle in the mesh is defined by three vertices (corner points) in 3D space.
Additionally, each triangle has a normal vector that indicates the direction of its surface, which helps determine the orientation of the triangle.


Users can also import STL files using Volmdlr and then use our features to perform powerful studies.

import STL files
***************

.. code-block:: python

    import volmdlr.cloud
    import volmdlr.core
    import volmdlr.stl as vmstl

    stl = vmstl.Stl.load_from_file('path/to/your/stl/file.stl')
    shell = stl.to_closed_shell()
    shell.alpha = 0.3
    shell.babylonjs()

If we break it down we have:

1. **Importing Required Modules:**

.. code-block:: python

    import volmdlr.cloud
    import volmdlr.core
    import volmdlr.stl as vmstl

In this section, the necessary modules from the volmdlr library are imported.
It seems to be organized into sub-modules such as cloud, core, and stl.
These modules likely provide different functionalities related to cloud computing, core geometry operations, and STL file manipulation, respectively.

2. **Loading a STL File:**

.. code-block:: python

    stl = vmstl.Stl.load_from_file('path/to/your/stl/file.stl')

Here, an STL file is loaded using the ``Stl.load_from_file()`` method from the vmstl sub-module.
The path to your STL file should be specified in the ``'path/to/your/stl/file.stl'`` placeholder.
The stl variable is assigned the resulting Stl object, which likely holds the geometry and other properties of the loaded STL file.

3. **Converting to ClosedShelll:**

.. code-block:: python

    shell = stl.to_closed_shell()

The ``to_closed_shell()`` method is called on the stl object.
This method presumably converts the loaded STL geometry into a closed shell object.
A closed shell is a common representation of a solid 3D object in computer graphics, where the outer surface completely encloses the interior volume.
The resulting shell variable holds this closed shell object.

5. **Generating Babylon.js Visualization:**

.. code-block:: python

    shell.babylonjs()

This line generates a 3D visualization of the closed shell using the Babylon.js format.
Babylon.js is a JavaScript framework for rendering 3D graphics in web browsers.

Simplify volume using Cloud of points
*************************************

Sometimes when working with stl files, it is needed to simplify the model to have a rough aproximation of its shape in order to run some space related calculations.
To achieve it, you can do as follows:

.. code-block:: python

    import volmdlr.cloud
    import volmdlr.core
    import volmdlr.stl as vmstl

    stl = vmstl.Stl.load_from_file('path/to/your/stl/file.stl')

    list_points = stl.extract_points_BIS()
    pointcloud3d = volmdlr.cloud.PointCloud3D(list_points)
    shell2 = pointcloud3d.to_shell(resolution=15)

    shell2.babylons()

In detail, each part of the code mean:

1. **Importing Required Modules and loading STL file:**

.. code-block:: python

    import volmdlr.cloud
    import volmdlr.core
    import volmdlr.stl as vmstl

    stl = vmstl.Stl.load_from_file('path/to/your/stl/file.stl')

As also shown in previous example, first we import all packages needed and then the STl file is loaded using the ``Vmstl.Stl.load_from_file`` method substituting ``'path/to/your/stl/file.stl'`` by the actual path to your stl file.

.. code-block:: python

    list_points = stl.extract_points_BIS()

The `extract_points_BIS()`` method is called on the stl object. This method extracts the points (vertices) from the STL geometry.
The resulting list_points variable holds the list of extracted points.

.. code-block:: python

    pointcloud3d = volmdlr.cloud.PointCloud3D(list_points)

A ``PointCloud3D`` object is created using the extracted list of points (list_points).
This step essentially converts the list of points into a point cloud representation suitable for further processing and visualization.

.. code-block:: python

    shell2 = pointcloud3d.to_shell(resolution=15)

The ``to_shell()`` method is called on the pointcloud3d object, converting the point cloud into a shell.
In this context, "shell" refers to a simplified surface representation of the point cloud. The resolution parameter determines how detailed the shell should be.
