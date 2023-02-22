from dessia_common.files import BinaryFile
from dessia_common.workflow.blocks import (CadView, ClassMethod,
                                           ClassMethodType, Export, MethodType,
                                           ModelMethod)
from dessia_common.workflow.core import Pipe, Workflow

import volmdlr as vm
from volmdlr import step

read_step_method_type = ClassMethodType(
    class_=vm.step.Step, name="from_stream"
)
cls_method_step = ClassMethod(read_step_method_type, name="Stepfile")

to_volumemodel_method_type = MethodType(
    class_=vm.step.Step, name="to_volume_model"
)
method_volumemodel = ModelMethod(
    to_volumemodel_method_type, name="VolumeModel"
)


cadview_block = CadView(name="Display3D")

export_html = Export(
    method_type=MethodType(vm.core.VolumeModel, "to_html_stream"),
    filename="VolumeModel.html",
    extension="html",
    text=True,
    name="Export_html",
)

pipes = [
    Pipe(cls_method_step.outputs[0], method_volumemodel.inputs[0]),
    Pipe(method_volumemodel.outputs[0], cadview_block.inputs[0]),
    Pipe(method_volumemodel.outputs[0], export_html.inputs[0]),
]


workflow = Workflow(
    [
        cls_method_step,
        method_volumemodel,
        cadview_block,
        export_html,
    ],
    pipes,
    method_volumemodel.outputs[0],
    name="From stl to volume model",
)

workflow.description = "Import step to Volume Model"

# workflow.plot()

Dir_file_component = "shell0.step"
with open(Dir_file_component, "rb") as stream_bin:
    bin_file = BinaryFile(Dir_file_component)
    bin_file.write(stream_bin.read())
    model = bin_file

input_values = {
    workflow.input_index(cls_method_step.inputs[0]): model,
}

workflow_run = workflow.run(input_values=input_values)

# workflow_run.output_value.babylonjs()
