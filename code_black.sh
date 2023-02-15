#!/bin/sh
if black --check -l120 .; then
    echo "Black found no formatting errors."
else
    echo "Black found formatting errors."
    echo "Please run 'black -l120 .' in source directory to format the code."
    echo "You can run 'black --check --diff --color -l120' the difference with expected formatting."
    exit 1
fi
