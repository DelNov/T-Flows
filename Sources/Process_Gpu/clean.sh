#!/bin/bash

find . -name "*.vtk" | xargs rm -f
find . -name "*.vtu" | xargs rm -f
