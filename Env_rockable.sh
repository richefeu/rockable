#!/bin/bash
# 
# Jean-Mathieu Vanson
# date: 10/08/2023
#
# Environment for Rockable exec

# Get the path of the install dir of rockable
pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}"
if ([ -h "${SCRIPT_PATH}" ]); then
  while([ -h "${SCRIPT_PATH}" ]); do cd `dirname "$SCRIPT_PATH"`;
  SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`pwd`;
popd  > /dev/null

# add rockable install dir to PATH
export PATH=${SCRIPT_PATH}/INSTALL:${PATH}
export LD_LIBRARY_PATH=${SCRIPT_PATH}/INSTALL:$LD_LIBRARY_PATH

