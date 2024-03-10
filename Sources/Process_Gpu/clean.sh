#!/bin/bash

find . -name "*.vtk" | xargs rm -f
find . -name "*.vtu" | xargs rm -f
find . -name "*.cfn" | xargs rm -f
find . -name "*.dim" | xargs rm -f
