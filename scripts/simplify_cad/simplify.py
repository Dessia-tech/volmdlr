import volmdlr.stl as vmstl
import volmdlr as vm
import volmdlr.simplify_cad
import math

stl_file = '../stl/a320_ENGINE_RIGHT.stl'
stl = vmstl.Stl.from_file(stl_file)
triangles = []
for k, tri in enumerate(stl.triangles):
    if not math.isclose(tri.area(), 0, abs_tol=1e-6):
        triangles.append(tri)

stl = vm.stl.Stl(triangles)
primitive = stl.to_closed_shell()

vol = vm.core.VolumeModel([primitive])
bb = vol.bounding_box
number = 20
size = max([(bb.xmax - bb.xmin) / number, (bb.ymax - bb.ymin) / number, (bb.zmax - bb.zmin) / number])
boxes_wrap = volmdlr.simplify_cad.BoxesWrap(vol, size)
# negative_boxe = boxes_wrap._generate_boxes_one_direction(vm.X3D, True)
all_boxes = boxes_wrap._generate_inside_boxes()
boxes = [boxes_wrap.define_block(i) for i in all_boxes]
primitive_merge = boxes_wrap._smart_glue_boxes(vm.Z3D, boxes)
primitive_merge.babylonjs()