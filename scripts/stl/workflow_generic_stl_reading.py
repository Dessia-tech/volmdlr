import dessia_common.workflow as wf
from dessia_common.typings import MethodType

import volmdlr as vm
from volmdlr import stl

read_stl_method_type = wf.blocks.ClassMethodType(class_=vm.stl.Stl, name="from_binary_stream")
cls_method_stl = wf.blocks.ClassMethod(read_stl_method_type, name="STLFile")

to_volumemodel_method_type = wf.blocks.MethodType(class_=vm.stl.Stl, name="to_volume_model")
method_volumemodel = wf.blocks.ModelMethod(to_volumemodel_method_type, name="VolumeModel")


cadview_block = wf.blocks.CadView(name="Display3D")

export_html = wf.blocks.Export(
    method_type=MethodType(vm.core.VolumeModel, "to_html_stream"),
    filename="VolumeModel.html",
    extension="html",
    text=True,
    name="Export_html",
)

pipes = [
    wf.Pipe(cls_method_stl.outputs[0], method_volumemodel.inputs[0]),
    # Display and export definition
    wf.Pipe(method_volumemodel.outputs[0], cadview_block.inputs[0]),
    wf.Pipe(method_volumemodel.outputs[0], export_html.inputs[0]),
]


workflow_stl = wf.core.Workflow(
    [
        cls_method_stl,
        method_volumemodel,
        cadview_block,
        export_html,
    ],
    pipes,
    method_volumemodel.outputs[0],
    name="From stl to volume model",
)

workflow_stl.description = "Import STL to Volume Model"

# workflow_stl.plot()


dict_workflow_stl = {i: j.name for i, j in enumerate(workflow_stl.inputs)}

# workflow_stl._check_platform()

# =============================================================================
# Platform insertion
# =============================================================================

# from dessia_api_client.users import PlatformUser
# platform = PlatformUser(api_url="https://api.platform.dessia.ovh")
# r = platform.objects.create_object_from_python_object(workflow_stl)
