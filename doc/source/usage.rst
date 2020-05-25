Usage
=====


VolumeModel
-----------
A model is set instanciating a VolumeModel:

.. autoclass:: volmdlr.VolumeModel

Prior to that, primitives have to be instanciated. See :ref:`primitives-label`.

FreeCAD binding
---------------

Once a VolumeModel is instanciated, one can call the FreeCADExport method:

.. autoclass:: volmdlr.VolumeModel
  :members: FreeCADExport

Sweep
-----
'Sweep' creates pipes with the following tool :

.. autoclass:: volmdlr.primitives3D.Sweep

To understand how Sweep works, see :ref:`primitives-label` in 3D advanced primitives.

Tutorials
---------

The scripts folder contains some examples of the capabilities of volmdlr:

https://github.com/Dessia-tech/volmdlr/tree/master/scripts
