import volmdlr as vm
import volmdlr.primitives3d as primitives3D


resolution = 0.0010

box = primitives3D.Block(
    vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.2, 0, 0),
                vm.Vector3D(0, 0.1, 0), vm.Vector3D(0, 0, 1)),
    alpha=0.6)

box_red = primitives3D.Block(
    vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.4, 0, 0),
                vm.Vector3D(0, 0.4, 0), vm.Vector3D(0, 0, 0.4)),
    color=(0.2, 1, 0.4), alpha=0.6)

box_red.color = (0.4, 0.1, 0.1)
box_red.name = 'box_red'

box_green = box_red.frame_mapping(vm.Frame3D(vm.Point3D(-0.4, 0, -0.1), vm.Vector3D(1, 0, 0),
                          vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'new', copy=True)

box_green.color = (0.1, 1, 0.1)
box_green.name = 'box_green'
shell1 = box_red.union(box_green)[0]
shell1.color = (0.2, 0.7, 0.47)
shell1.alpha = 0.6
vm.core.VolumeModel([shell1]).babylonjs()
#
box_blue = box.frame_mapping(vm.Frame3D(vm.Point3D(0.1, 0, 0), vm.Vector3D(1, 0, 0),
                          vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'old', copy=True)
box_blue.color = (0.1, 0.1, 1)
box_blue.name = 'box_blue'

box_blue2 = box.frame_mapping(vm.Frame3D(vm.Point3D(0.3, 0, 0), vm.Vector3D(1, 0, 0),
                          vm.Vector3D(0, 1.8, 0), vm.Vector3D(0, 0, 1)), 'old', copy=True)
box_blue2.color = (0.1, 0.1, 1)
box_blue2.name = 'box_blue2'
shell2 = box_blue.union(box_blue2)[0]
shell2.color = (0.5, 0.5, 0.5)
shell2.alpha = 0.6
vm.core.VolumeModel([shell2]).babylonjs()
union_box = shell1.union(shell2)
subtraction_box = shell1.subtract(shell2)
intersection_box = shell1.intersection(shell2)
subtraction_closedbox = shell1.subtract_to_closed_shell(shell2)


for new_box in [union_box, subtraction_box, subtraction_closedbox, intersection_box]:
    for shell in new_box:
        shell.color = (1, 0.1, 0.1)
        shell.alpha = 0.6
    vm.core.VolumeModel(new_box + [shell1, shell2]).babylonjs()


contour = vm.wires.ClosedPolygon2D([vm.Point2D(0, 0), vm.Point2D(-1, 0),
                                    vm.Point2D(-1, 1), vm.Point2D(1, 1),
                                    vm.Point2D(1, -1), vm.Point2D(0, -1)])
extrude1 = vm.primitives3d.ExtrudedProfile(vm.O3D, vm.X3D, vm.Y3D,
                                           contour, [],  -vm.Z3D)
extrude1.color = (0.4, 0.1, 0.1)
extrude1.alpha = 0.6
frame1 = vm.Frame3D(vm.Point3D(0, 0, 0.5), 2*vm.X3D, 2*vm.Y3D, vm.Z3D)
block1 = vm.primitives3d.Block(frame1)
block1.color = (0.2, 0.7, 0.47)
block1.alpha = 0.6
union1 = extrude1.union(block1)
union1[0].alpha = 0.6
union1[0].color = (1, 0.1, 0.1)
vol = vm.core.VolumeModel([union1[0], block1, extrude1])
vol.babylonjs()

extrude2 = vm.primitives3d.ExtrudedProfile(vm.O3D, vm.X3D, vm.Y3D,
                                           contour, [], -vm.Z3D)
extrude2.color = (0.4, 0.1, 0.1)
extrude2.alpha = 0.6
frame1 = vm.Frame3D(vm.Point3D(0, 0.5, -0.5), 2*vm.X3D, vm.Y3D, vm.Z3D)
block1 = vm.primitives3d.Block(frame1)
block1.color = (0.2, 0.7, 0.47)
block1.alpha = 0.6

bbox1 = extrude2.bounding_box
bbox2 = block1.bounding_box
union1 = extrude2.union(block1)
union1[0].alpha = 0.6
union1[0].color = (1, 0.1, 0.1)
vol = vm.core.VolumeModel([union1[0]])
vol.babylonjs()

frame1 = vm.Frame3D(vm.Point3D(0, 0, 0), vm.X3D, vm.Y3D, vm.Z3D)
frame2 = vm.Frame3D(vm.Point3D(0, 0.8, 0), vm.X3D, vm.Y3D, vm.Z3D)
block1 = vm.primitives3d.Block(frame1)
block2 = vm.primitives3d.Block(frame2)
block3 = block1.union(block2)
block3[0].alpha = 0.6
block3[0].color = (1, 0.1, 0.1)
vol = vm.core.VolumeModel([block1, block2, block3[0]])
vol.babylonjs()

frame2 = vm.Frame3D(vm.Point3D(0, 1, 0), vm.X3D, vm.Y3D, 3*vm.Z3D)
block1 = vm.primitives3d.Block(frame1)
block2 = vm.primitives3d.Block(frame2)
block3 = block1.union(block2)
block3[0].alpha = 0.6
block3[0].color = (1, 0.1, 0.1)
vol = vm.core.VolumeModel([block1, block2, block3[0]])
vol.babylonjs()