import os, glob
import time
import nbformat
from nbconvert import PythonExporter


# Converting .ipynb to .py
def convert_notebook(notebook_path, module_path):
    with open(notebook_path) as fh:
        nb = nbformat.reads(fh.read(), nbformat.NO_CONVERT)

    exporter = PythonExporter()
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True
    source, meta = exporter.from_notebook_node(nb)
    lines = source.splitlines()
    with open(module_path, 'w+') as fh:
        for line in lines:
            # Removing the prompt commands from notebook
            if line[:20] != 'get_ipython().system' and '.plot(' not in line and line[:1] != '!'  and line[:12] != 'from IPython' and line[:4] != 'HTML':
                fh.write(line + '\n')


tutorials = []

# Converting all tutorials automatically
for tutorial in glob.glob("*.ipynb"):
    tutorials.append(tutorial[:-6])

for tutorial in tutorials:
    convert_notebook(tutorial + '.ipynb', tutorial + '.py')

for i in range(len(tutorials)):
    tutorials[i] += '.py'

# Testing if all tutorials exists before launching them
for tutorial_name in tutorials:
    if not os.path.isfile(tutorial_name):
        raise FileNotFoundError(f'Script {tutorial_name} does not exists in CI tutorials')

# Executing tutorials
print('Executing tutorials for CI:')
total_time = time.time()
top_level_dir = os.getcwd()
times = {}
for tutorial_name in tutorials:
    print(f'\t* {tutorial_name}')
    # Reset dir
    os.chdir(top_level_dir)
    # Change cwd
    if '/' in tutorial_name:
        tutorial_folder = '/'.join(tutorial_name.split('/')[:-1])
        if tutorial_folder:
            tutorial_folder = os.path.join(top_level_dir, tutorial_folder)
            os.chdir(tutorial_folder)
    file_name = tutorial_name.split('/')[-1]
    t = time.time()
    with open(file_name, 'r', encoding='utf-8') as tutorial:
        exec(tutorial.read())
    t = time.time() - t
    times[tutorial_name] = t

print('Computation times:')
for tutorial_name, t in sorted(times.items(), key=lambda x: x[1]):
    print(f'* tutorial {tutorial_name}: {round(t, 3)} seconds ')
    
total_time = time.time() - total_time
print(f'Total time for CI tutorials: {total_time}')
