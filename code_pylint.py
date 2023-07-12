import os
import random
import sys
from datetime import date, timedelta
import math

from pylint import __version__
from pylint.lint import Run

MIN_NOTE = 9.6

UNWATCHED_ERRORS = ['fixme', 'trailing-whitespace', 'import-error', 'missing-final-newline']

EFFECTIVE_DATE = date(2023, 1, 31)

WEEKLY_DECREASE = 0.03

MAX_ERROR_BY_TYPE = {
                     "wrong-spelling-in-comment": 199,
                     "wrong-spelling-in-docstring": 143,
                     'invalid-name': 184,
                     'no-member': 3,
                     'inconsistent-return-statements': 4,
                     'unused-variable': 22,
                     'arguments-differ': 62,
                     'too-many-locals': 75,
                     'unused-argument': 32,
                     'too-many-arguments': 34,
                     'line-too-long': 12,
                     'too-many-branches': 26,
                     'too-many-statements': 18,
                     'super-init-not-called': 15,
                     'no-name-in-module': 11,
                     'abstract-method': 32,
                     'duplicate-code': 10,
                     'arguments-renamed': 52,
                     'too-many-ancestors': 24,
                     'too-few-public-methods': 2,
                     'too-many-public-methods': 12,
                     'too-many-instance-attributes': 15,
                     'protected-access': 4,
                     'undefined-loop-variable': 2,
                     'unspecified-encoding': 1,
                     'too-many-function-args': 4,
                     'too-many-nested-blocks': 7,
                     'too-many-return-statements': 3,
                     'cyclic-import': 1,
                     'undefined-variable': 0,  # 2 when gmsh is fixed
                     'broad-except': 1,
                     "broad-exception-caught": 2,
                     'too-many-boolean-expressions': 2,
                     'too-many-lines': 4,
                     'consider-using-with': 1,
                     'unnecessary-dunder-call': 2,
                     'chained-comparison': 2,
                     'consider-using-generator': 1,
                     'import-outside-toplevel': 5,
                     'unsubscriptable-object': 1,
                     'signature-differs': 1}

ERRORS_WITHOUT_TIME_DECREASE = ["too-many-locals", "too-many-branches", "too-many-arguments", "too-many-statements",
                                "too-many-nested-blocks", "too-many-instance-attributes", "abstract-method",
                                "no-name-in-module", "too-many-public-methods", "too-many-ancestors",
                                "protected-access", "cyclic-import", "line-too-long", "too-many-lines", "no-member",
                                "too-few-public-methods", "duplicate-code", "too-many-return-statements",
                                "import-outside-toplevel", "arguments-differ", "arguments-renamed",
                                "too-many-boolean-expressions"]

limit_time_effect = False
if os.environ.get('DRONE_BRANCH', '') in ['master', 'testing']:
    limit_time_effect = True
    print(f"Limiting time effect of 21 days as we are on {os.environ['DRONE_BRANCH']}")

if os.environ.get('DRONE_TARGET_BRANCH', '') in ['master', 'testing']:
    limit_time_effect = True
    print(f"Limiting time effect of 21 days as we are targetting {os.environ['DRONE_TARGET_BRANCH']}")

if limit_time_effect:
    EFFECTIVE_DATE += timedelta(days=21)


print("pylint version: ", __version__)

time_decrease_coeff = 1 - (date.today() - EFFECTIVE_DATE).days / 7.0 * WEEKLY_DECREASE

f = open(os.devnull, "w")

old_stdout = sys.stdout
sys.stdout = f

results = Run(["volmdlr", "--output-format=json", "--reports=no"], do_exit=False)
# `exit` is deprecated, use `do_exit` instead
sys.stdout = old_stdout

PYLINT_OBJECTS = True
if hasattr(results.linter.stats, "global_note"):
    pylint_note = results.linter.stats.global_note
    PYLINT_OBJECT_STATS = True
else:
    pylint_note = results.linter.stats["global_note"]
    PYLINT_OBJECT_STATS = False


def extract_messages_by_type(type_):
    return [m for m in results.linter.reporter.messages if m.symbol == type_]


error_detected = False
error_over_ratchet_limit = False

if PYLINT_OBJECT_STATS:
    stats_by_msg = results.linter.stats.by_msg
else:
    stats_by_msg = results.linter.stats["by_msg"]

print(f'Errors / Allowed errors: {sum(stats_by_msg.values())} / {sum(MAX_ERROR_BY_TYPE.values())})')

for error_type, number_errors in stats_by_msg.items():
    if error_type not in UNWATCHED_ERRORS:
        base_errors = MAX_ERROR_BY_TYPE.get(error_type, 0)

        if error_type in ERRORS_WITHOUT_TIME_DECREASE:
            max_errors = base_errors
        else:
            max_errors = math.ceil(base_errors * time_decrease_coeff)

        time_decrease_effect = base_errors - max_errors
        # print('time_decrease_effect', time_decrease_effect)

        if number_errors > max_errors:
            error_detected = True
            print(
                f"\nFix some {error_type} errors: {number_errors}/{max_errors} "
                f"(time effect: {time_decrease_effect} errors)")

            messages = extract_messages_by_type(error_type)
            messages_to_show = sorted(random.sample(messages, min(30, len(messages))), key=lambda m: (m.path, m.line))
            for message in messages_to_show:
                print(f"{message.path} line {message.line}: {message.msg}")
        elif number_errors < max_errors:
            print(f"\nYou can lower number of {error_type} to {number_errors+time_decrease_effect}"
                  f" (actual {base_errors})")

for error_type in MAX_ERROR_BY_TYPE:
    if error_type not in stats_by_msg:
        print(f"You can delete {error_type} entry from MAX_ERROR_BY_TYPE dict")

if error_detected:
    raise RuntimeError("Too many errors\nRun pylint volmdlr to get the errors")

if error_over_ratchet_limit:
    raise RuntimeError("Please lower the error limits in code_pylint.py MAX_ERROR_BY_TYPE according to warnings above")

print("Pylint note: ", pylint_note)
if pylint_note < MIN_NOTE:
    raise ValueError(f"Pylint not is too low: {pylint_note}, expected {MIN_NOTE}")

print("You can increase MIN_NOTE in pylint to {} (actual: {})".format(pylint_note, MIN_NOTE))
