
import volmdlr as vm
import volmdlr.faces
import volmdlr.primitives3d

primitive1 = vm.primitives3d.ExtrudedProfile.load_from_file('test1.json')
primitive2 = vm.primitives3d.ExtrudedProfile.load_from_file('test2.json')
vol = vm.core.VolumeModel([primitive2, primitive1])
vol.babylonjs()

open_shell = primitive1.subtract(primitive2)
open_shell[0].babylonjs()


