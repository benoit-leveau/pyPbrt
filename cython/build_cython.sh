#!/bin/bash

# compile
python setup.py build_ext --build-lib=build/pypbrt --build-temp=.build --cython-c-in-temp

# install
for folder in build/pypbrt/* ; do touch $folder/__init__.py ; done
touch build/pypbrt/__init__.py

