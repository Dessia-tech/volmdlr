import pydocstyle
import os
from glob import glob
import random
from datetime import date

file_list = filter(lambda z: not z.endswith("__init__.py"),
                   [y for x in os.walk('./volmdlr')
                    for y in glob(os.path.join(x[0], '*.py'))])

UNWATCHED_ERRORS = [
    # Do not watch these errors
    'D100', 'D104', 'D105', 'D107',
    'D200', 'D202', 'D203', 'D204', 'D206', 'D210', 'D212',
    'D301', 'D302',
    'D401', 'D402', 'D407', 'D408', 'D409',
    'D412', 'D415', 'D418'
]

MAX_ERROR_BY_TYPE = {
    # http://www.pydocstyle.org/en/stable/error_codes.html
    'D100': 1,
    'D101': 67,
    'D102': 713,
    'D103': 35,
    'D104': 1,
    'D105': 1,
    'D106': 1,
    'D107': 1,

    'D200': 1,
    'D201': 1,
    'D202': 1,
    'D203': 1,
    'D204': 1,
    'D205': 279,
    'D206': 1,
    'D207': 1,
    'D208': 7,
    'D209': 1,
    'D210': 1,
    'D211': 1,
    'D212': 1,
    'D213': 2,
    'D214': 1,
    'D215': 1,

    'D300': 6,
    'D301': 1,
    'D302': 1,

    'D400': 377,
    'D401': 1,
    'D402': 1,
    'D403': 85,
    'D404': 6,
    'D405': 1,
    'D406': 1,
    'D407': 1,
    'D408': 1,
    'D409': 1,
    'D410': 1,
    'D411': 1,
    'D412': 1,
    'D413': 3,
    'D414': 1,
    'D415': 1,
    'D416': 1,
    'D417': 8,
    'D418': 1,
}

error_detected = False
error_over_ratchet_limit = False
ratchet_limit = 9
effective_date = date(2022, 11, 28)
today = date.today()
weekly_decrease = 5
time_decrease = int((today - effective_date).days/7. * weekly_decrease)


code_to_errors = {}
for error in pydocstyle.check(file_list, ignore=UNWATCHED_ERRORS):
    code_to_errors.setdefault(error.code, [])
    code_to_errors[error.code].append(error)

code_to_number = {code: len(errors) for code, errors in code_to_errors.items()}

for error_code, number_errors in code_to_number.items():
    if error_code not in UNWATCHED_ERRORS:
        max_errors = max(MAX_ERROR_BY_TYPE.get(error_code, 0) - time_decrease, 0)

        if number_errors > max_errors:
            error_detected = True
            print(f'\nFix some {error_code} errors: {number_errors}/{max_errors}')

            errors = code_to_errors[error_code]
            errors_to_show = sorted(random.sample(errors, min(30, len(errors))),
                                    key=lambda m: (m.filename, m.line))
            for error in errors_to_show:
                print(f'{error.filename} line {error.line}: {error.message}')
        elif max_errors - ratchet_limit <= number_errors < max_errors:
            print(f'\nYou can lower number of {error_code} to {number_errors + time_decrease} (actual {max_errors + time_decrease})')
        elif number_errors < max_errors - ratchet_limit:
            error_over_ratchet_limit = True
            print(f'\nYou MUST lower number of {error_code} to {number_errors + time_decrease} (actual {max_errors + time_decrease})')

if error_detected:
    raise RuntimeError('Too many errors\nRun pydocstyle volmdlr to get the errors')

if error_over_ratchet_limit:
    raise RuntimeError('Please lower the error limits in code_pydocstyle.py MAX_ERROR_BY_TYPE according to warnings above')
