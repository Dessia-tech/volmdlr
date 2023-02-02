#!/bin/bash
# check pep8 formatting for changed files

PEP8_CMD_TO_RUN='python3 -m autopep8 -i $(git diff --cached --name-only --diff-filter MARC | grep -E \.py$)'

CHANGED_FILES=$(git diff --cached --name-only --diff-filter MARC | grep -E \.py$)
if [[ -z "$CHANGED_FILES" ]]
  then exit 0
fi

DETECTED_CHANGES=$(python3 -m autopep8 -d $CHANGED_FILES)
if [[ -n "$DETECTED_CHANGES" ]]
  then
  echo -e "\npep8 non conforming changes detected, please run :\n"
  echo -e "\t$PEP8_CMD_TO_RUN\n"
  echo -e "& stage your changes"
  exit 1
fi

exit 0
