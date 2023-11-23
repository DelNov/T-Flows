---
project: Generate
output_dir: ./Html
author: T-Flows Development Team
email: "your.email@example.com"
preprocessor: gfortran -cpp -E
project_github: https://github.com/DelNov/T-Flows
project_download: https://github.com/DelNov/T-Flows/archive/refs/heads/main.zip
src_dir: ./Temporary
---

<!-- This is a comment -->
<!-- <div style="text-align: center;">                                                        -->
<!-- <p><img alt="T-Flows Logo" src="../../../logo_100_percent.png" title="Program Logo"></p> -->
<!-- </div>                                                                                   -->

## Introduction

The sub-program _Generate_ is a utility for generating a computational grid
to be used in numerical simulations with T-Flows.  The grid generation process
includes creating structured or block-structured meshes and performing
unstructured cell refinements where necessary. _Generate_ also handles domain
connectivity, periodicity, and potentially smoothing regions within the mesh.

### Features

_Generate_ was developed in the 1990s by one of the members of the development
team due to inavailability of mesh generation softwares (we are talking about
the age when two-dimensional simulations were a standard and three-dimensional
simulations were a luxury, no open source mesh generators were available for
three-dimensional grids and internet itself was in its infancy). Today,
_Generate_ is considered outdated and is kept only for compatibility with older
simulation cases. T-Flows users are strongly encouraged to use
[GMSH](https://gmsh.info) or any other free or commercial mesh generator
capable to export ANSYS' ```.msh``` or ```.cas``` format (or even the legacy
GAMBIT's ```.neu```) file format and then use _Convert_ sub-program to convert
the grids to T-Flows' format.

Despite being obsolete, _Generate_ is useful for maintaining continuity with
legacy simulation cases, ensuring that results can still be reproduced or
compared with newer cases.

### Getting Started

...

