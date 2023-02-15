#!/bin/sh
# check black formatting on all files

if black --check -l120 .; then
  echo ""
  echo "Black found no formatting errors."
else
  echo ""
  echo "Black found formatting errors."
  echo "Please run 'black -l120 .' in source directory to format the code."
  echo "You can run 'black --check --diff --color -l120' to see the difference with expected formatting."
  echo "Consider upgrading black to latest version if you're still getting an error after formatting the code."
  exit 1
fi
