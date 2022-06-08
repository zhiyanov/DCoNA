#!/bin/bash

mkdir build
cd build
cmake ../dcona/native
make
cp *.so ../dcona/core/extern/
cd ..
rm -rf build
