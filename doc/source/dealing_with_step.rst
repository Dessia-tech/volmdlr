====
STEP
====

A STEP file, which stands for "Standard for the Exchange of Product model data," is a standardized file format used for exchanging 3D CAD (Computer-Aided Design) models and related data between different software applications and systems.
It is a widely used format in the field of engineering, manufacturing, and product design.

A STEP file is designed to facilitate the seamless exchange of complex 3D product information across various software platforms.
It provides a standardized way to represent geometric and non-geometric data associated with a product's design, such as its shape, dimensions, assembly structure, materials, and other attributes.

To help users coming from another CAD software, the Volmdlr library provides a module called `step.py that enables you to import STEP files and perform powerfull analysis within volmdlr features.

Import a STEP file
*****************

Here's a example of how you can do it:


.. code-block:: python

    from volmdlr import step

    stepfile = step.Step.from_file('/path/to/your/step/file.step')

    volume = stepfile.to_volume_model()
    volume.primitives[0].color = (1, 0.2, 0.2)
    volume.primitives[0].alpha = 0.6
    volume.babylonjs()

If we break it down we have:

1. **Importing Required Modules:**

.. code-block:: python

    from volmdlr import step

The script imports the step module from the volmdlr library.
This module provides functionality to work with STEP (Standard for the Exchange of Product model data) files, which are used for exchanging 3D CAD models and associated information.

2. **Loading a STEP File:**

.. code-block:: python

    stepfile = step.Step.from_file('/path/to/your/step/file.step')

This line loads a STEP file using the ``Step.from_file()`` method.
You need to replace ``'/path/to/your/step/file.step'`` with the actual file path to the STEP file you want to load.
This file contains the 3D CAD model and its data.

3. **Converting to Volume Model:**

.. code-block:: python

    volume = stepfile.to_volume_model()

The loaded STEP file is converted into a volume model using the ``to_volume_model()`` method.
A volume model represents the geometry of the CAD model in volmdlr.
This conversion allows you to manipulate and visualize the 3D geometry.

4. **Adjusting Visual Properties:**

.. code-block:: python

    volume.primitives[0].color = (1, 0.2, 0.2)
    volume.primitives[0].alpha = 0.6

This section modifies the visual properties of the first primitive in the volume model.
The color property is set to (1, 0.2, 0.2), which represents a reddish hue using RGB values.
The alpha property is set to 0.6, introducing a level of transparency to the object.

5. **Generating Babylon.js Visualization:**

.. code-block:: python

    volume.babylonjs()

This line generates a 3D visualization of the volume model using the Babylon.js format.
Babylon.js is a JavaScript framework for rendering 3D graphics in web browsers.
This step prepares the data and structure needed to render the 3D object, considering the visual properties and
modifications applied earlier.

In summary, this script loads a 3D CAD model from a STEP file, converts it into a volume model, modifies the visual
appearance of the model's first primitive (color and transparency), and then generates a 3D visualization using the
Babylon.js format. The resulting visualization displays the modified CAD model with the specified color and transparency settings.
