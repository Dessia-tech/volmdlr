# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- WIP : read STEP files

### Changed
- MPLPlot now returns only ax instead of fig, ax 

## [v0.1.11]

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
