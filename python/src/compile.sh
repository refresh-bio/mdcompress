#!/bin/bash

cd ../nanobind
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DNANOBIND_BUILD_TESTS=OFF \
         -DNANOBIND_BUILD_EXAMPLES=OFF \
         -DNANOBIND_BUILD_STUBS=OFF
make -j
cd ../../src

rm *.so
g++ -shared -fPIC -O3 -std=c++20 python_bind.cpp \
    -I ../../src/lib_mdc \
    -L ../../bin \
    -I ../nanobind/include \
    $(python3-config --includes) \
    -L ../nanobind/build/tests \
    -o mdc$(python3-config --extension-suffix) \
    -lnanobind-static -Wl,--whole-archive -lmdc -Wl,--no-whole-archive
python3 example.py

