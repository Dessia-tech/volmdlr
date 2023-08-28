==================
Usefull operations
==================

Export to json file
******************

You can export any volmdlr objects to a json file. To do so, you can just call the method ``save_to_file``.

Here is an example how you can save a plan face 3d to a .json file:

.. code-block:: python

    import volmdlr
    from volmdlr import surfaces, faces

    face = faces.PlaneFace3D.from_surface_rectangular_cut(surfaces.Plane3D(volmdlr.OXYZ), -1, 1, -1, 1)
    face.save_to_file('path/to/where/to/save/your/file.json')


Import from json file
********************

If you have a json file containing volmdlr objects, you can use the method ``dessia_common.core.DessiaObect.load_from_file`` to import it.

Example:

.. code-block:: python

    import dessia_common

    volmdlr_object = dessia_common.core.DessiaObject.load_from_file('path/to/your/file.json')

Matplotlib plots
****************

To have a matplotlib visualization of a volmdlr, you can call the plot method in any object. See the following example
#tobe continued
Model Visualization
*******************

To Have a 3D visulaization of your model, You can use the babylonjs() method in any 3D object.
#tobe continued
