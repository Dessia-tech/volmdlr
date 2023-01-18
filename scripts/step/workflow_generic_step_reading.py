import io

import dessia_common.workflow as wf
from dessia_common.typings import MethodType

import volmdlr as vm
from volmdlr import step

read_step_method_type = wf.blocks.ClassMethodType(
    class_=vm.step.Step, name="from_stream"
)
cls_method_step = wf.blocks.ClassMethod(read_step_method_type, name="Stepfile")

to_volumemodel_method_type = wf.blocks.MethodType(
    class_=vm.step.Step, name="to_volume_model"
)
method_volumemodel = wf.blocks.ModelMethod(
    to_volumemodel_method_type, name="VolumeModel"
)


cadview_block = wf.blocks.CadView(name="Display3D")

export_html = wf.blocks.Export(
    method_type=MethodType(vm.core.VolumeModel, "to_html_stream"),
    filename="VolumeModel.html",
    extension="html",
    text=True,
    name="Export_html",
)

pipes = [
    wf.Pipe(cls_method_step.outputs[0], method_volumemodel.inputs[0]),
    wf.Pipe(method_volumemodel.outputs[0], cadview_block.inputs[0]),
    wf.Pipe(method_volumemodel.outputs[0], export_html.inputs[0]),
]


workflow = wf.core.Workflow(
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

workflow.plot()


model = io.FileIO("shell0.step", "r")

input_values = {
    workflow.input_index(cls_method_step.inputs[0]): model,
}

workflow_run = workflow.run(input_values=input_values)

# workflow_run.output_value.babylonjs()
