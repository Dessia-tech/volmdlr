#!/bin/sh
# check black formatting on all files

if black --check .; then
  echo ""
  echo "Black found no formatting errors."
else
  echo ""
  echo "Black found formatting errors."
  echo "Please run 'black .' in source directory to format the code."
  echo "You can run 'black --check --diff --color .' to see the difference with expected formatting."
  echo ""
  echo "Consider upgrading black to latest version if you're still getting an error after formatting the code."
  echo ""
  echo "Note: Black is configurated in 'pyproject.toml' to have a 120 character line length."
  echo "It is equivalent than using the 'l120' argument: 'black -l120 .'"
  exit 1
fi
