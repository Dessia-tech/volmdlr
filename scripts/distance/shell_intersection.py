import volmdlr as vm
import volmdlr.primitives3d as primitives3d

resolution = 0.0010

box = primitives3d.Block(
    vm.Frame3D(vm.Point3D(0.1, 0.1, 0.1), vm.Vector3D(0.15, 0.0, 0.0),
               vm.Vector3D(0.0, 0.15, 0.0), vm.Vector3D(0.0, 0.0, 0.15)),
    alpha=0.6)

box.frame_mapping_inplace(vm.Frame3D(vm.Point3D(-0.07, -0.07, -0.07),
                             vm.Point3D(1, 0, 0),
                             vm.Point3D(0, 1, 0),
                             vm.Point3D(0, 0, 1)), side='old')

box_red = primitives3d.Block(
    vm.Frame3D(vm.Point3D(-0.04, -0.04, -0.04), vm.Vector3D(0.1, 0.0, 0.0),
               vm.Vector3D(0.0, 0.1, 0.0), vm.Vector3D(0.0, 0.0, 0.1)),
    color=(0.2, 1, 0.4), alpha=0.6)

assert type(box_red.shell_intersection(box, resolution=0.001)) == tuple

vol = vm.core.VolumeModel([box, box_red])
vol.babylonjs()
