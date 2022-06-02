#!/bin/bash

mkdir build
cd build
cmake ../native
make
cp *.so ../core/extern/
cd ..
rm -rf build
