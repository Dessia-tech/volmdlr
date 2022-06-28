#!/bin/bash
# check pep8 formatting for all files

PEP8_CMD_TO_RUN='python3 -m autopep8 -i volmdlr/*.py'

DETECTED_CHANGES=$(python3 -m autopep8 -d volmdlr/*.py)
if [[ -n "$DETECTED_CHANGES" ]]
  then
  echo -e "\npep8 non conforming changes detected, please run :\n"
  echo -e "\t$PEP8_CMD_TO_RUN\n"
  echo -e "& stage your changes before pushing"
  exit 1
fi

exit 0
