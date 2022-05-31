import os
import time

unittests = [
    # faces
    'faces/test_closedshell3d.py',
    'faces/test_planeface3d.py'
]

# Testing if all unittests exists before launching them
for script_name in unittests:
    if not os.path.isfile(script_name):
        raise FileNotFoundError(f'Script {script_name} does not exists in CI unittests')

# Executing unittests
print('Executing unittests for CI:')
total_time = time.time()
top_level_dir = os.getcwd()
times = {}
for script_name in unittests:
    print(f'\t* {script_name}')
    # Reset dir
    os.chdir(top_level_dir)
    # Change cwd
    if '/' in script_name:
        test_folder = '/'.join(script_name.split('/')[:-1])
        if test_folder:
            script_folder = os.path.join(top_level_dir, test_folder)
            os.chdir(script_folder)
    file_name = script_name.split('/')[-1]
    t = time.time()
    with open(file_name, 'r', encoding='utf-8') as unittest:
        exec(unittest.read())
    t = time.time() - t
    times[script_name] = t

print('Computation times:')
for script_name, t in sorted(times.items(), key=lambda x: x[1]):
    print(f'* script {script_name}: {round(t, 3)} seconds ')

total_time = time.time() - total_time
print(f'Total time for CI unittests: {total_time}')
