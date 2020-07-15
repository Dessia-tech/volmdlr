Tutorial
=====

This guide can help you using tools developed with the Volume Modeler Volmdlr.

Models
------
In volmdlr, you can create Models as a face contained in a plane, Cylinder (CylindricalFace3D),
a Tore (ToroidalFace3D) and more (BSpline, ConicalFace, SphericalFace WIP).
Please refers to :ref:`primitives-label` in 3D advanced primitives to know how it works.

To show the final result, we use Babylonjs.

VolumeModel
-----------
A model is set instanciating a VolumeModel:

.. autoclass:: volmdlr.VolumeModel

Prior to that, primitives have to be instanciated. See :ref:`primitives-label`.

FreeCAD binding (Optional)
---------------

Once a VolumeModel is instanciated, one can call the FreeCADExport method:

.. autoclass:: volmdlr.VolumeModel
  :members: FreeCADExport

Sweep
-----
'Sweep' creates pipes with the following tool :

.. autoclass:: volmdlr.primitives3D.Sweep

To understand how Sweep works, see :ref:`primitives3D-label` in 3D advanced primitives.

Shell3D
-------

If you want to create a simple Face3D or more, you can create them thanks to our tools see :ref:`primitives3D-label` in 3D advanced primitives.

:Example: 
  >>> You can create a Shell3D : 
	import volmdlr as vm
	(create faces with our examples of each Face3D)
	shell = vm.Shell3D([face1, face2, etc..])
	volumodel = vm.VolumeModel([shell])
	volumodel.babylonjs() (to show it in babylonjs interface)

More Tutorials
--------------

The scripts folder contains some examples of the capabilities of volmdlr:

https://github.com/Dessia-tech/volmdlr/tree/master/scripts 