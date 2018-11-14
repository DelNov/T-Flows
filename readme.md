#T-Flows 

T-Flows (stands for Turbulent Flows) is a Computational Fluid Dynamics (CFD) program featuring second order accurate finite volume discretization of incompressible Navier-Stokes equations with heat transfer and species transport.  It is written in Fortran 90 and uses Message Passing Interface (MPI) for parallel execution.

At the moment a few Reyolds Stress Navier-Stokes (RANS) models are implemented:
- Classical k-eps, 
- Elliptically blended k-eps-zeta-f, 
- Reynolds Stress Model (RSM) from Manceau and Hanjalic and 
- RSM from Jakirlic and Hanjalic.

Each model can be ran in unstready mode, thus effectivelly becoming Unsteady RANS (URANS) model.

In additiont to the RANS models outlined above, T-Flows also features Large Eddy Simulation (LES) with following Sub-Grid Scale (SGS) models:
- Smagorinsky,
- Germano's dynamic model,
- Wall-Adapting Local Eddy viscosity (WALE) model.

The code also features a hybrid RANS/LES model, based on k-eps-zeta-f close to the wall and Germano's dynamic model far from the wall.  

Code has four sub-programs at the moment and uses own proprietary format for computational grids (.geo, .cns).
- Generate - to generate mesh from ASCII .dom to own format.
- Convert - to convert mesh from .neu, .cgns to own format.
- Divide - to partition the mesh for parallel execution
- Process - to solve discretized Navier-Stokes equations.

Build depends on the following tools: 
ld, make, gfortran, git and mpi (only for parallel execution).

Build tested with following tools and versions:
- ld >= 2.24
- gfortran >= 4.8
- mpich >= 3.0.4
- openmpi >= 1.6.5
- make >= 3.81
- git >= 1.91

To download current version of the code use:

git clone https://github.com/DelNov/T-Flows/

Compiled programs are in Binaries/ folder.
Sources are in Sources/Generate/, Sources/Convert/, Sources/Divide/, Sources/Process/ .

Tests/ folder contains a few cases for RANS, LES and hybrid methods.
