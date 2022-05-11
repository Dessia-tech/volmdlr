import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.primitives3d as p3d

sphere_surface = vmf.SphericalSurface3D(vm.Frame3D(vm.O3D, vm.X3D, vm.Y3D, vm.Z3D), 0.3)
sphere_surface.plot()

face = sphere_surface.rectangular_cut(0, 2.5, 0, 0.5)
face.babylonjs()
face.plot()

# Try to generate a sphere from revoltion
radius = 0.03
p1 = vm.Point2D(0, -radius)
p2 = vm.Point2D(radius, 0.)
p3 = vm.Point2D(0, radius)
arc = vme.Arc2D(p1, p2, p3)
line = vme.LineSegment2D(p3, p1)
contour = vmw.Contour2D([arc, line])
sphere_revolution = p3d.RevolvedProfile(vm.X3D, vm.X3D, vm.Y3D, contour,
                                        vm.X3D, vm.Y3D, angle=1.3)
sphere_revolution.babylonjs()


sphere = p3d.Sphere(vm.O3D, 0.06)
sphere.babylonjs()