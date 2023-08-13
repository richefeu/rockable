#!/bin/bash
#
# created on 07/08/2023
# Author: JM Vanson

echo $*
my_dir=$PWD

CMAKE_OPTION=-DCMAKE_INSTALL_PREFIX=${my_dir}/INSTALL

if [[ "$*" == *"-pbc"* ]] 
then
   CMAKE_OPTION=${CMAKE_OPTION}" -DROCKABLE_ENABLE_PERIODIC=ON"
fi

if [[ "$*" == *"-clean"* ]] 
then
   rm -r ${my_dir}/BUILD
   rm -r ${my_dir}/INSTALL
fi

mkdir BUILD INSTALL
cd BUILD
cmake ${CMAKE_OPTION} ..
make -j
make install
