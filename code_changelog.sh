#!/bin/bash

lines=$(git diff origin/"$DRONE_TARGET_BRANCH"..origin/"$DRONE_SOURCE_BRANCH" -- CHANGELOG.md README.md CONTRIBUTING.md --unified=0| wc -l)
echo "$lines lines modified on CHANGELOG.md or other doc files in PR $DRONE_SOURCE_BRANCH -> $DRONE_TARGET_BRANCH"


if [ "$lines" -eq 0 ]
  then
  echo -e "\nCHANGELOG.md has not been updated. Update it for the PR to be accepted in CI.\n"
  exit 1
fi
