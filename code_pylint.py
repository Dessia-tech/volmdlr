import os
import sys

import random
from pylint import __version__
from pylint.lint import Run

MIN_NOTE = 9.17

UNWATCHED_ERRORS = ['fixme', 'trailing-whitespace', 'import-error', 'missing-final-newline']


MIN_NOTE = 7.9

MAX_ERROR_BY_TYPE = {
                     'invalid-name': 1023,
                     'no-else-return': 74,
                     'consider-using-f-string': 155,
                     'no-member': 127,
                     'inconsistent-return-statements': 18,
                     'unused-variable': 61,
                     'arguments-differ': 44,
                     'too-many-locals': 75,
                     'line-too-long': 32,
                     'unused-argument': 48,
                     'too-many-arguments': 65,
                     'line-too-long errors': 32,
                     'consider-using-enumerate': 22,
                     'too-many-branches': 27,
                     'too-many-statements': 21,
                     'super-init-not-called': 19,
                     'no-name-in-module': 22,
                     'abstract-method': 24,
                     'empty-docstring': 18,
                     'duplicate-code': 11,
                     'no-self-use': 16,
                     'arguments-renamed': 17,
                     'too-many-ancestors': 17,
                     'expression-not-assigned': 12,
                     'non-parent-init-called': 11,
                     'too-few-public-methods': 10,
                     'too-many-public-methods': 8,
                     'use-implicit-booleaness-not-comparison': 8,
                     'too-many-instance-attributes': 12,
                     'protected-access': 8,
                     'undefined-loop-variable': 5,
                     'unspecified-encoding': 9,
                     'too-many-function-args': 8,
                     'too-many-nested-blocks': 6,
                     'attribute-defined-outside-init': 9,
                     'too-many-return-statements': 4,
                     'consider-merging-isinstance': 6,
                     'cyclic-import': 4,
                     'consider-iterating-dictionary': 2,
                     'raise-missing-from': 4,
                     'no-else-raise': 4,
                     'no-else-continue': 4,
                     'undefined-variable': 4,
                     'no-else-break': 4,
                     'unnecessary-list-index-lookup': 4,
                     'simplifiable-if-expression': 3,
                     'redefined-builtin': 3,
                     'broad-except': 3,
                     'too-many-boolean-expressions': 3,
                     'too-many-lines': 3,
                     'redundant-keyword-arg': 3,
                     'no-value-for-parameter': 1,
                     'c-extension-no-member': 2,
                     'access-member-before-definition': 2,
                     'modified-iterating-list': 2,
                     'consider-using-with': 2,
                     'consider-using-get': 2,
                     'unnecessary-dunder-call': 2,
                     'unnecessary-lambda': 2,
                     'eval-used': 2,
                     'chained-comparison': 2,
                     'missing-module-docstring': 2,
                     'unbalanced-tuple-unpacking': 2,
                     'bad-staticmethod-argument': 1,
                     'consider-using-generator': 1,
                     'use-maxsplit-arg': 1,
                     'wildcard-import': 1,
                     'cell-var-from-loop': 1,
                     'import-outside-toplevel': 1,
                     'unsubscriptable-object': 1,
                     # No tolerance errors
                     'unidiomatic-typecheck': 0,
                     'unexpected-special-method-signature': 0,
                     'bare-except': 0,
                     'function-redefined': 0,
                     'superfluous-parens': 0,
                     'unnecessary-comprehension': 0,
                     'unused-wildcard-import': 0,
                     'wrong-import-order': 0,
                     'ungrouped-imports': 0,
                     'assignment-from-none': 0,
                     'non-iterator-returned': 0,
                     'consider-using-max-builtin': 0,
                     'consider-using-min-builtin': 0,
                     'format-string-without-interpolation': 0,
                     'bad-indentation': 0,
                     'implicit-str-concat': 0,
                     'return-in-init': 0,
                     'redefined-outer-name': 0,
                     'no-self-argument': 0,
                     'consider-using-from-import': 0,
                     'duplicate-string-formatting-argument': 0,
                     'wrong-import-position': 0,
                     'singleton-comparison': 0,
                     'unreachable': 0,
                     'consider-using-in': 0,
                     'unused-import': 0
                     }

print('pylint version: ', __version__)

f = open(os.devnull, 'w')

old_stdout = sys.stdout
sys.stdout = f

results = Run(['volmdlr', '--output-format=json', '--reports=no'], do_exit=False)
# `exit` is deprecated, use `do_exit` instead
sys.stdout = old_stdout

PYLINT_OBJECTS = True
if hasattr(results.linter.stats, 'global_note'):
    pylint_note = results.linter.stats.global_note
    PYLINT_OBJECT_STATS = True
else:
    pylint_note = results.linter.stats['global_note']
    PYLINT_OBJECT_STATS = False


def extract_messages_by_type(type_):
    return [m for m in results.linter.reporter.messages if m.symbol == type_]


error_detected = False
error_over_ratchet_limit = False

if PYLINT_OBJECT_STATS:
    stats_by_msg = results.linter.stats.by_msg
else:
    stats_by_msg = results.linter.stats['by_msg']

for error_type, number_errors in stats_by_msg.items():
    if error_type not in UNWATCHED_ERRORS:
        max_errors = MAX_ERROR_BY_TYPE.get(error_type, 0)

        if number_errors > max_errors:
            error_detected = True
            print(f'\nFix some {error_type} errors: {number_errors}/{max_errors}')

            messages = extract_messages_by_type(error_type)
            messages_to_show = sorted(random.sample(messages, min(30, len(messages))),
                                      key=lambda m: (m.path, m.line))
            for message in messages_to_show:
                print(f'{message.path} line {message.line}: {message.msg}')
        elif number_errors < max_errors:
            print(f'\nYou can lower number of {error_type} to {number_errors} (actual {max_errors})')


if error_detected:
    raise RuntimeError('Too many errors\nRun pylint volmdlr to get the errors')

if error_over_ratchet_limit:
    raise RuntimeError('Please lower the error limits in code_pylint.py MAX_ERROR_BY_TYPE according to warnings above')

print('Pylint note: ', pylint_note)
# if pylint_note > MIN_NOTE + RATCHET_NOTE:
#     raise ValueError(f'MIN_NOTE in code_pylint.py is too low, increase to at least {MIN_NOTE + RATCHET_NOTE}, max {pylint_note}')
if pylint_note < MIN_NOTE:
    raise ValueError(f'Pylint not is too low: {pylint_note}, expected {MIN_NOTE}')

print('You can increase MIN_NOTE in pylint to {} (actual: {})'.format(pylint_note,
                                                                      MIN_NOTE))
