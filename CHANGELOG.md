# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unrealeased
* Duplicated methods cleaned


*improvements to union operations

### Fixed
* Remove Copy param from movement of primitives and add inplace methods



* improvements to union operations

### Fixed
Remove Copy param from movement of primitives and add inplace methods
* Wire3D babylonjs


## v0.4.0
### Fixed
* various fixes in cuts of wires and contours
* Fix of missing face in Union
* following dessia_common v0.7.0


## v0.3.0

### New Features
* Bspline with dimensions
* cut_by_line for Surface2D
* Bspline merge

### Fixed
* Various Steps improvement
* Bspline periodicity in step reading
* sewing improvements
* Substraction of shells

## v0.2.10

### New Features

* union of shells (only with planeface for the moment 
* Sewing of polygon3D
* Concav hull of PointCloud2D

## v0.2.9

### New Features

* support STL import & export
* point cloud2D & cloud3D

## v0.2.8

### New Features

* support stringIO in step save

### Fixes

* depack of point2D
* to_vector2D

### Performance improvements

* better bounding box for cylindrical face


## [v0.2.7]
### Changed
* direction vector of linesegments are now normalized

### New Features

* straight line area for BsplineCurve2D
* split of circleby start end
* closedpolygon2d is_trigo
* Auto-adaptative camera/edge width babylonjs
* splitting of bsplinecurve2d
* BezierSurface3D implemented
* added rotation and translation for faces
* new classes BezierCurve2D and BezierCurve3D
* spherical surface
* (core): update plot_data method
* update plot_data methods in wires and edges
* step almost working for cylindrical, conical toroidal
* difference between intersections and crossings
* plot_data version set to 0.3.8 or above

### Fixes

* support of mixed vector point in to step
* remove debug mode babylonjs
* remove sci notation in step export
* use stable cdn for babylonjs
* sweep extrusion length
* line circle intersection with tolerance, normal and dir vector for arc
* offset of wire
* remove useless non serializable attr
* secondmoment area from straight lines
* reversed faces in extrusion correction
* enhancement of rotation/translation of shells
* bug fix BezierCurve2D and 3D
* eq and hash for basis and frames
* shell and frame mapped shell correctly read
* small try except added for step reading
* all SHAPE_REPRESENTATION are now read
* Arc3D from step full debug
* arc3d to 2d in bspline3d surface
* missing faces at end of sweep
* splitting faces and arcs
* perf in display nodes and toroidal aspect
* setup.py requires plot_data>=0.3.9
* (primitives2d): serialization
* debug of shell method
* porting shells methods
* Debug of conical faces
* Porting cylinders and hollow
* porting from missing from_contour3d for planeface
* reading steps, but artefact on faces
* Correcting arc from_step

### Performance improvements

* LineSegment2D.points is non serializable attribute
* ClosedPolygon2D.line_segment is non_serializable_attributes
* Optimization of mesh generation

#### Refactorings
* (edges): put data argument back into Arc2D.plot_data()
* (edges): redefined Arc2D.plot_data()

## v0.2.6

### Changed
- debugs on frame 2D 

### Optimized
- babylon data generation speed up

## v0.2.5

### Added
- translation and rotation for various primitives

### Changed
- Frame3D rotation takes also into account origin
- following plot_data v0.5.3

## v0.2.4
### Added
- handle spherical surfaces
- positionning of parts in STEP reading

## v0.2.1
### Added
- step export

## v0.2

### Changed
- modules *2D or *3D renamed in *2d, *3d
- point and vector declared with their x, y, z vm.Point2D((0, 0)) -> vm.Point2D(0, 0)
- separating in new modules: display, wires, edges...
- PEP8: method names
- PointAtCurvilinearAbscissa changed to point_at_abscissa
- MPLPlot changed to plot()
- plot now returns only ax instead of fig, ax 

## v0.1.11

### Added 
- Calculate the distance between LineSegment3D/LS3D, Arc3D/LS3D, Arc3D/Arc3D and between CylindricalFace3D too.
- Use PlaneFace3D with contours2D in a classic way and use it with contours3D with a 'from_contours3d' as CylindricalFace3D does.
- Calculate the distance between CylindricalFace3D and PlaneFace3D.
- Calculate the distance between CylindricalFace3D, PlaneFace3D and ToroidalFace3D.
- contours2d.tessel_points which gives all points of a contour2d, and .points the end points of primitives.
- Implementation of ConicalFace3D in Core and RevolvedProfile.
- Implementation of SphericalFace3D in Core.
- BSplineFace3D works.

### Changed
- cut_contours in Face3D which take all points from a Contour2D, not one side like before. Furthermore, it is light and quick.

## [v0.1.10]
- typings
- workflow to instanciate point

## [v0.1.9]

### Added
- mesh module

## [v0.1.8]

### Added
- color and alpha options for various primitives
- line segments intersection
 
### Debug
- arcs: is_trigo and angle were sometimes false

## [v0.1.7]

### Added
- random vector and points
- dashed line option in babylon of LineSegment3D
- Measure2D
- babylon_data: a dict language to describe models to be unpacked by a babylonjs unpacker

### Removed
- constants o2D, x2D, y2D...: use O2D, X2D...

### Changed
- Mesure -> Measure3D
