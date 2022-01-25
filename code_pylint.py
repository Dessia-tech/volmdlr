from pylint.lint import Run

MIN_NOTE = 5.03

MAX_ERROR_BY_TYPE = {'cyclic-import': 6,
                     'too-many-lines': 3,
                     'bare-except': 0,
                     'no-else-return': 184,
                     'no-self-use': 17,
                     'no-member': 185,
                     'unexpected-special-method-signature': 0,
                     'too-many-locals': 89,
                     'too-many-nested-blocks': 7,
                     'inconsistent-return-statements': 108,
                     'arguments-differ': 90,
                     'too-many-arguments': 64,
                     'undefined-variable': 3,
                     'function-redefined': 0,
                     'attribute-defined-outside-init': 9,
                     'simplifiable-if-expression': 3,
                     'redefined-builtin': 3,
                     'unnecessary-comprehension': 0,
                     'consider-using-enumerate': 30,
                     'no-value-for-parameter': 2,
                     'abstract-method': 7,
                     'wildcard-import': 1,
                     'unused-wildcard-import': 0,
                     'too-many-return-statements': 6,
                     'eval-used': 2,
                     'too-many-statements': 28,
                     'superfluous-parens': 0,
                     'chained-comparison': 2,
                     'wrong-import-order': 8,
                     'unused-variable': 95,
                     'unused-import': 9,
                     'super-init-not-called': 28,
                     'consider-using-f-string': 173,
                     'too-many-branches': 30,
                     'consider-merging-isinstance': 6,
                     'too-many-instance-attributes': 10,
                     'unused-argument': 79,
                     'undefined-loop-variable': 10,
                     'consider-using-with': 3,
                     'use-maxsplit-arg': 1,
                     'broad-except': 3,
                     'consider-iterating-dictionary': 4,
                     'raise-missing-from': 4,
                     'unspecified-encoding': 9,
                     'import-outside-toplevel': 3,
                     'consider-using-get': 2,
                     'ungrouped-imports': 0,
                     'bad-staticmethod-argument': 1,
                     'arguments-renamed': 14,
                     'too-many-public-methods': 12,
                     'consider-using-generator': 1
                     }

import os
import sys
f = open(os.devnull, 'w')

old_stdout = sys.stdout
sys.stdout = f

results = Run(['volmdlr', '--output-format=json', '--reports=no'], do_exit=False)
# `exit` is deprecated, use `do_exit` instead
sys.stdout = old_stdout

pylint_note = results.linter.stats.global_note
print('Pylint note: ', pylint_note)
assert pylint_note >= MIN_NOTE
print('You can increase MIN_NOTE in pylint to {} (actual: {})'.format(pylint_note,
                                                                      MIN_NOTE))


def extract_messages_by_type(type_):
    return [m for m in results.linter.reporter.messages if m.symbol == type_]


uncontrolled_errors = {}
error_detected = False
for error_type, number_errors in results.linter.stats.by_msg.items():
    if error_type in MAX_ERROR_BY_TYPE:
        if number_errors > MAX_ERROR_BY_TYPE[error_type]:
            error_detected = True
            print('Fix some {} errors: {}/{}'.format(error_type,
                                                     number_errors,
                                                     MAX_ERROR_BY_TYPE[error_type]))
            for message in extract_messages_by_type(error_type):
                print('{} line {}: {}'.format(message.path, message.line, message.msg))
        elif number_errors < MAX_ERROR_BY_TYPE[error_type]:
            print('You can lower number of {} to {} (actual {})'.format(
                error_type, number_errors, MAX_ERROR_BY_TYPE[error_type]))

    else:
        if not error_type in uncontrolled_errors:
            uncontrolled_errors[error_type] = number_errors

if uncontrolled_errors:
    print('Uncontrolled errors', uncontrolled_errors)

if error_detected:
    raise RuntimeError('Too many errors\nRun pylint dessia_common to get the errors')
