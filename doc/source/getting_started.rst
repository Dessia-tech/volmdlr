Volmdlr: the absolute basics for beginners
------------------------------------------
Welcome to the absolute beginner’s guide to volmdlr! If you have comments or suggestions, please don’t hesitate to reach out!

Welcome to volmdlr!
^^^^^^^^^^^^^^^^^^^

The volmdlr library is an open-source Python library primarily developed by Dessia Technologies,
aimed at providing 3D modeling capabilities based on Boundary Representation (B-Rep) purely in Python.
The library is designed to be easy to use, efficient, and customizable.

Boundary Representation is a technique used in computer-aided design (CAD) and computer-aided
manufacturing (CAM) systems to represent the geometry of a solid object as a collection of surfaces
and curves. B-Rep is widely used in CAD and CAM due to its ability to accurately represent complex
shapes and to provide a rich set of operations for manipulating and analyzing these shapes.

Installing volmdlr through pip (Windows and Linux users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This package requires Python 3.9 or above. Please follow the instructions
below to install the package. Depending on your needs, you can choose between
two types of installation.

To install the latest version of the package you need to run the following
command::

  pip install volmdlr
  # or
  pip3 install volmdlr

To install a specific version of the package you would issue the following
command::

  pip install volmdlr==0.1.0
  # or
  pip3 install volmdlr==0.1.0

Developer installation
^^^^^^^^^^^^^^^^^^^^^^

First, clone the package. Then, enter the newly created volmdlr repository. Finally, develop the setup.py file, and you are good to go ! ::

  git clone https://github.com/Dessia-tech/volmdlr.git

  cd volmdlr

  python3 setup.py develop --user
  # or whatever version you are using :
  python3.x setup.py develop --user

Requirements
^^^^^^^^^^^^

The installation of volmdlr requires the installation of other packages listed
in the file setup.py and in the table below. These libraries will be
automatically installed when you install volmdlr.

=============  ===============  ===========
Dependency     Minimum Version  Usage
=============  ===============  ===========
packaging          latest       computation
dessia_common      >=0.10.0     computation      
Cython             latest       computation
numpy              latest       computation
matplotlib         latest       display
scipy              latest       computation
geomdl             latest       computation
networkx           latest       computation
triangle           latest       computation
plot_data          >=0.10.9     display
kaitaistruct       latest       computation
binaryornot        latest       computation
sympy              latest       computation
trimesh            latest       computation
rtree              latest       computation
gmsh               latest       computation
=============  ===============  ===========

Troubleshooting
^^^^^^^^^^^^^^^

If the installation is successful but your IDE don't recognize the package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case you may have several versions of Python installed on your
computer. Make sure the `pip` command points to the right Python version, or
that you have selected the desired Python version in your IDE.
You can force the installation of the package on a given Python version by
executing this command::

  python -m pip install volmdlr

You have to specify the Python version you are working with by replacing
`python` by the Python of your choice. For example, `python3`, `python3.8`,
`python3.9`, etc.
