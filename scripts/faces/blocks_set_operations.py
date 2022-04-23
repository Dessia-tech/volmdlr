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
                          vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'new')

box_green.color = (0.1, 1, 0.1)
box_green.name = 'box_green'
shell1 = box_red.union(box_green)[0]
shell1.merge_faces()
shell1.color = (0.2, 0.7, 0.47)
shell1.alpha = 0.6
vm.core.VolumeModel([shell1]).babylonjs()
#
box_blue = box.frame_mapping(vm.Frame3D(vm.Point3D(0.1, 0, 0), vm.Vector3D(1, 0, 0),
                          vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'old')
box_blue.color = (0.1, 0.1, 1)
box_blue.name = 'box_blue'

box_blue2 = box.frame_mapping(vm.Frame3D(vm.Point3D(0.3, 0, 0), vm.Vector3D(1, 0, 0),
                          vm.Vector3D(0, 1.8, 0), vm.Vector3D(0, 0, 1)), 'old')
box_blue2.color = (0.1, 0.1, 1)
box_blue2.name = 'box_blue2'
shell2 = box_blue.union(box_blue2)[0]
shell2.merge_faces()
shell2.color = (0.5, 0.5, 0.5)
shell2.alpha = 0.6
vm.core.VolumeModel([shell2]).babylonjs()
union_box = shell1.union(shell2)
union_box[0].merge_faces()
subtraction_box = shell1.subtract(shell2)
intersection_box = shell1.intersection(shell2)
intersection_box[0].merge_faces()
subtraction_closedbox = shell1.subtract_to_closed_shell(shell2)
subtraction_closedbox[0].merge_faces()


for new_box in [union_box, subtraction_box, subtraction_closedbox, intersection_box]:
    for shell in new_box:
        shell.color = (1, 0.1, 0.1)
        shell.alpha = 0.6
    vm.core.VolumeModel(new_box + [shell1, shell2]).babylonjs()


block1 = primitives3D.Block(frame=vm.Frame3D(origin=vm.Point3D(0.2, 0.6, 0.2),
                                                  u=vm.Point3D(0.4, 0, 0),
                                                  v=vm.Point3D(0, 1.2, 0),
                                                  w=vm.Point3D(0, 0, 0.4)),
                            color=(1, 0.2, 0.1), alpha=0.6)
block2 = primitives3D.Block(frame=vm.Frame3D(origin=vm.Point3D(0.8+0.4, 0.45 + 0.15 , 0.1),
                                                  u=vm.Point3D(1.6, 0, 0),
                                                  v=vm.Point3D(0, 0.9, 0),
                                                  w=vm.Point3D(0, 0, 0.2)),
                            color=(0.2, 0.2, 0.2), alpha=0.6)

block3 = primitives3D.Block(frame=vm.Frame3D(origin=vm.Point3D(0.8+0.4, 0.45 + 0.15 , 0.05),
                                                  u=vm.Point3D(1.6, 0, 0),
                                                  v=vm.Point3D(0, 0.9, 0),
                                                  w=vm.Point3D(0, 0, 0.2)),
                            color=(0.2, 0.2, 0.2), alpha=0.6)

for battery_layout_max_enclosing_box in [block1.union(block2)[0], block1.union(block3)[0]]:
    battery_layout_max_enclosing_box.merge_faces()
    battery_layout_max_enclosing_box.color = (1, 0.2, 0.2)
    battery_layout_max_enclosing_box.alpha = 0.6
    battery_layout_max_enclosing_box.babylonjs()