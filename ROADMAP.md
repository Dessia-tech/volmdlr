# Volmdlr Road map

This is the change roadmap for volmdlr.
Once implemented, update CHANGELOG.md

## Short term (next month)
### Bugs 
* Edge width too big when rendering far from the origin

### Naming
Formalize methods names: 
* discretise, discretization_points, polygon_points

### Refactorings
* rotation, translation and frame_mapping: Divide each in two methods,
where one modifies the object inplace and the other creates a new object. Delete copy parameter
* replace normalize() methodes by to_normal(). The object should no more be modified inplace.
* edges: CompositePrimitive Inheritance (use UML graphs)
* duplicated methods
* Grid of point class

### New features
 * BSpline3D extrusion

## Midterm (6 months)
### New features
* ClosedShell3D offset
* Remove triangle dependency: Surface2d with inner contour triangulation uses triangle.
Wheel unavailable for windows for python3.9. 

## LongTerm (1 year)