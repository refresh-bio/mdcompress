#!/bin/bash

#mkokot_TODO: change this such that there is no need to specify so path with PYTHONPATH=...

#this is for VS Code: -DCMAKE_EXPORT_COMPILE_COMMANDS=ON


rm -rf build
#cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

cmake --build build -- -j$(nproc)
PYTHONPATH=build/ python3 src/example.py

