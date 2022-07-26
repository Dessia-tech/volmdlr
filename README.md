<h1 align="center">
  <img src="https://partage.dessia.tech/thumbnail/7861a783126742be8fe8/1024/Logo_Dessia_transparent_web.png" style="width:300px"><br/>Volmdlr
</h1>

<h4 align="center">
  A computations-oriented python VOLume MoDeLeR with STEP support for import and export
</h4>

<div align="center">
  <a href="http://dessia.tech/"><img src="https://img.shields.io/website-up-down-green-red/http/dessia.tech.svg"></a>  
  <a href="https://GitHub.com/Dessia-tech/volmdlr/stargazers/"><img src="https://badgen.net/github/stars/Dessia-tech/volmdlr"></a>  
  <a href="https://drone-opensource.dessia.tech/Dessia-tech/volmdlr"><img src="https://drone-opensource.dessia.tech/api/badges/Dessia-tech/volmdlr/status.svg?branch=master"></a>
  <a href="https://pypi.org/project/volmdlr/"><img src="https://img.shields.io/pypi/v/volmdlr.svg"></a>
  <a href="https://github.com/Dessia-tech/volmdlr/graphs/contributors"><img src="https://img.shields.io/github/contributors/Dessia-tech/volmdlr.svg"></a>
  <a href="https://github.com/Dessia-tech/volmdlr/issues"><img src="https://img.shields.io/github/issues/Dessia-tech/volmdlr.svg"></a>
</div>

<div align="center">
  <a href="#description"><b>Description</b></a> |
  <a href="#features"><b>Features</b></a> |
  <a href="#user-installation"><b>User Installation</b></a> |
  <a href="#dev-installation"><b>Dev Installation</b></a> |
  <a href="https://github.com/Dessia-tech/volmdlr/tree/master/scripts"><b>Usage</b></a> |
  <a href="https://documentation.dessia.tech/volmdlr/"><b>Documentation</b></a> |
  <a href="#licence"><b>Licence</b></a> |
  <a href="#contributors"><b>Contributors</b></a> |
</div>

## Description

Volmdlr is a python volume modeler used as a CAD platform.
It is simple to understand and operate.
With it, you can easily create 3D models.
Check the exemples to see what you can do with this library.

<p align="center"><img src="https://raw.githubusercontent.com/Dessia-tech/volmdlr/master/doc/source/images/casing.jpg" width="40%" /> <img src="https://raw.githubusercontent.com/Dessia-tech/volmdlr/master/doc/source/images/casing_contours.png" width="55%" /></p>
<i>A casing is defined by a 2D contour formed with the primitive RoundedLineSegment2D. This contour is offset by the casing width.</i><br/><br/><br/>

<p align="center"><img src="https://raw.githubusercontent.com/Dessia-tech/volmdlr/master/doc/source/images/sweep1.jpg" width="45%" /> <img src="https://raw.githubusercontent.com/Dessia-tech/volmdlr/master/doc/source/images/sweepMPLPlot.jpg" width="50%" /></p>
<i>A Sweep is pipes, created with Circle2D/Arc2D which is contained in a Contour2D. You have to create the neutral fiber, i.e., the pipe’s road, with the primitive RoundedLineSegment3D.</i><br/><br/><br/>

<p align="center"><img src="https://raw.githubusercontent.com/Dessia-tech/volmdlr/master/doc/source/images/polygon.jpg" width="47%" /></p>
<i>A polygon is defined out of points. Random points are sampled and the tested whether they are inside or outside the polygon. They are plotted with the Matplotlib binding MPLPlot with custom styles:
- red if they are outside,
- blue if they are inside
</i><br/><br/><br/>

<p align="center"><img src="https://raw.githubusercontent.com/Dessia-tech/volmdlr/master/doc/source/images/bspline_surface_split.png" width="47%" /></p>
<i>A 3D B-spline surface splitted by a 3D B-spline curve.</i><br/><br/><br/>

## Features

- [x] Generate 2D and 3D geometries from python
- [x] Handles complexe geometries : B-spline curves and surfaces
- [x] Primitives provide computational tasks : distances, belonging, union, intersections, etc.
- [x] STEP/STL imports and exports
- [x] Geometries display in your web browser with [babylon.js](https://www.babylonjs.com/)

## User Installation

```bash
pip install volmdlr
# or
pip3 install volmdlr
```

## Dev Installation

Before using Volmdlr, be sure to have a C/C++ compiler (not necessary on Linux).  
N.B : With Windows you have to download one and allow it to read Python’s code.

First, [clone](https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories) the package.
Then, enter the newly created volmdlr repository.
Finally, develop the setup.py file, and you are good to go !

```bash
git clone https://github.com/Dessia-tech/volmdlr.git

cd volmdlr

python3 setup.py develop --user
# or whatever version you are using :
python3.x setup.py develop --user
```

## Usage

See the [script](https://github.com/Dessia-tech/volmdlr/tree/master/scripts) folder for examples

## Documentation

https://documentation.dessia.tech/volmdlr/

## Licence

100% opensource on LGPL licence. See LICENCE for more details.

## Contributors

- [DessiA team](https://github.com/orgs/Dessia-tech/people)
