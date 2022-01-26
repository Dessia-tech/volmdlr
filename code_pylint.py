from pylint.lint import Run

MIN_NOTE = 5.03

MAX_ERROR_BY_TYPE = {
                     'no-member': 185,
                     'no-else-return': 178,
                     'consider-using-f-string': 167,
                     'inconsistent-return-statements': 108,
                     'unused-variable': 95,
                     'arguments-differ': 90,
                     'too-many-locals': 89,
                     'unused-argument': 72,
                     'too-many-arguments': 64,
                     'consider-using-enumerate': 30,
                     'too-many-branches': 30,
                     'too-many-statements': 28,
                     'super-init-not-called': 28,
                     'abstract-method': 18,
                     'duplicate-code': 16,
                     'no-self-use': 16,
                     'arguments-renamed': 14,
                     'too-many-public-methods': 12,
                     'too-many-instance-attributes': 10,
                     'undefined-loop-variable': 10,
                     'unused-import': 9,
                     'unspecified-encoding': 9,
                     'too-many-nested-blocks': 9,
                     'attribute-defined-outside-init': 7,
                     'too-many-return-statements': 6,
                     'consider-merging-isinstance': 6,
                     'cyclic-import': 6,
                     'consider-iterating-dictionary': 4,
                     'raise-missing-from': 4,
                     'undefined-variable': 3,
                     'simplifiable-if-expression': 3,
                     'redefined-builtin': 3,
                     'broad-except': 3,
                     'too-many-lines': 3,
                     'redundant-keyword-arg': 3,
                     'no-value-for-parameter': 2,
                     'consider-using-with': 2,
                     'consider-using-get': 2,
                     'eval-used': 2,
                     'chained-comparison': 2,
                     'unbalanced-tuple-unpacking': 2,
                     'bad-staticmethod-argument': 1,
                     'consider-using-generator': 1,
                     'use-maxsplit-arg': 1,
                     'wildcard-import': 1,
                     'cell-var-from-loop': 1,
                     # No tolerance errors
                     'import-outside-toplevel': 0,
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
                     'consider-using-in': 0
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
