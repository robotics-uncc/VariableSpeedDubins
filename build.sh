#!/bin/bash

# if it does not exist, create a build folder
mkdir -p build
cd build

# run cmake
cmake ..

# run make and return to root directory
make
cd ..

