# T-Flows 

T-Flows (stands for Turbulent Flows) is a Computational Fluid Dynamics (CFD) program, originally developed at the Delft University of Technology, The Netherlands,  featuring second order accurate unstructured finite volume discretization of incompressible Navier-Stokes equations with heat transfer and species transport.  It is written in Fortran 90 and uses Message Passing Interface (MPI) for parallel execution.

Features: 
•	Unstructured (cell-centered) grid, arbitrary cell shapes 2nd-order accuracy; 
•	Linear solvers: preconditioned Conjugate Cradient (CG), Bi-Conjugate-Gradient (BiCG) and Conjugate Gradient Squared (CGS)
•	Basic and advanced RANS models
•	Smagorinsky and dynamic subgrid-scale models for LES
•	Hybrid RANS/LES methods
•	Conventional and advanced treatment of wall boundary conditions (wall integration, wall functions with roughness and compound wall treatment - blending wall function and wall integration approach)
•	Compatible with several open-source and commercial grid generators

At the moment the following models of turbulence are implemented:
•	Linear Eddy-Viscosity Models k-ε (standard high-Re, low-Re version of Abe, Kondoh and Nagano (AKN) with compound wall treatment)
•	Elliptic-Relaxation Eddy-viscosity Model - linear ζ − f with compound wall treatment
•	Second-moment Closures Full Re-stress models
-  Elliptic Blending model, 
-  Hanjalic-Jakirlic model
-  Hybrid Linear Eddy-Viscosity Models - Second-moment Closures model (Basara-Jakirlic)
•	LES (Large-eddy Simulations)
 - Standard and Dynamic Smagorinsky model
 - Wall-Adapting Local Eddy viscosity (WALE) model
•	Hybrid LES/RANS
 - DES
 - hybrid LES/RANS ζ − f model

All RANS models can be ran in unsteady mode, thus effectively becoming Unsteady RANS (URANS) model.

Code has four sub-programs at the moment and uses own proprietary format for computational grids (.geo, .cns).
- Generate -  use to generate mesh from ASCII .dom to own format.
- Convert  -  use to convert mesh from .neu, .cgns to own format.
- Divide   -  use to partition the mesh for parallel execution.
- Process  -  use to solve discretized Navier-Stokes equations.

Compilation of the main program and sub-programs depends on the following tools: 
ld, make, gfortran, git and mpi (only for parallel execution).

Compilation was tested with the following tools and versions:
- ld >= 2.24
- gfortran >= 4.8
- mpich >= 3.0.4
- openmpi >= 1.6.5
- make >= 3.81
- git >= 1.91

To download the current version of the code use:

 git clone https://github.com/DelNov/T-Flows/ `

Compiled programs are placed in the Binaries/ folder.

Sources are located in
Sources/Generate/
Sources/Convert/
Sources/Divide/
Sources/Process/

Tests/ folder contains files needed to run several standard benchmark cases for RANS, LES and hybrid methods. Xmgrace sources with referent DNS or experimental data are provided for the most of the test cases, as well as user functions for post-processing of data.

The program produces files in vtu format that can be open by ParaView open-source software for postprocessing. 
