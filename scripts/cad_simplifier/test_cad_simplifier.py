import os
import volmdlr.step as STEP
import volmdlr.core as vmc
from volmdlr import cad_simplification


# files_source_folder = '../step'
files_source_folder = ''

file_names = [
]


def print_seprator():
    print('============================================================================')


def simpflier_form():
    print_seprator()
    print('which simplifier method would you like to use? (type corresponding number)')
    print('1 - TrippleExtrusionSimplify')
    print('2 - OctreeBlockSimplify')
    chosen_option = int(input('Choose simplifier: '))
    print_seprator()
    return chosen_option


for filename in file_names:
    model = STEP.Step.from_file(filepath=os.path.join(files_source_folder, filename))
    volume_model = model.to_volume_model()
    simplifier_option = simpflier_form()
    if simplifier_option == 1:
        simplified_model = cad_simplification.TrippleExtrusionSimplify(volume_model).simplify()
    elif simplifier_option == 2:
        precision = int(input('How precise do you want the simplificatio to be (0-6)? (0  being not precise \n'
                              'at all and 6 closer to original)'))
        print_seprator()
        simplified_model = cad_simplification.OctreeBlockSimplify(volume_model.primitives[0]).simplify(precision)
    else:
        raise ValueError('Not valid choice. Please select one of the available options.')
    simplified_model.alpha = 0.6
    simplified_model.color = (1, .1, 0.1)
    both_model = vmc.VolumeModel(volume_model.primitives + [simplified_model])
    both_model.babylonjs()
