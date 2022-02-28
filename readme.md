# T-Flows 
<br />
1. [Introduction](#intro)
2. [Software requirements](#soft_req)
    1. [Minimim](#soft_req_min)
    2. [Higlhy desirable](#soft_req_des)
    3. [Optional](#soft_req_opt)
3. [Obtaining the code](#obtaining)
4. [Compiling the code](#compiling)
<br />


# Introduction <a name="intro"></a>
<br />
T-Flows is a computational fluid dynamics (CFD) program for simulation of turbulent, single and multiphase flows.  Numerical method is based on collocated finite volume method on unstructured arbitrary grids and turbulence models include a range of Reynolds-averaged Navier-Stokes (RANS) models, large eddy simulations (LES), as well as hybrid RANS-LES approach.  A more comprehensive list of turbulence models can be found [here](./Documentation/Manual/turbulence_models.md).

Multiphase models include an algebraic volume of fluid (VOF) method and Lagrangian particle tracking model.  Three-phase flows situations (two fluid phases with VOF and one solid phase as particles) are also supported.

> **_Note:_** In T-Flows, the Navier-Stokes equations are discretized in their _incompressible_ form, meaning _only_ that pressure and temperatures are _not_ linked through an equation of state.  All physical properties in T-Flows can be variable, but you should keep in mind that variable density does not mean compressibility.
<br />


# Software requirements <a name="soft_req"></a>
<br />

## Minimum software requirements <a name="soft_req_min"></a>
<br />

The bare minimum to get T-Flows running entails:

- make utility
- Fortran 2008 compiler
- standard C compiler

T-Flows is almost entirely written in Fortran 2008 (only one function is written in C) and the compilation is controlled by makefiles.  So, the the requirements listed above are a bare minimum for you to start using the code.  
t
Although there is, in principle, no restriction on the operating system on which you can use T-Flows, its natural habitat is Linux, as we develop test it on Linux, and Linux meets the minimum software requirements either _out of the box_, or with minimum installation effort.

> **_Note:_** We do not specify the minimum version for any of the required or recommended software.  We believe that if you are reading these pages, you do have access to a relatively recent hardware which also implies an up to date operating system and the associated tools.  
<br />

## Highly desirable software requirements <a name="soft_req_des"></a>
<br />

Although without meeting the minimum software requirements listed above you will not get anywhere, they alone will not get you very far either.  To make a practical use of T-Flows, it is highly desirable that you also have the following:

- [GMSH](https://gmsh.info)
- any other free or commercial mesh generator exporting ANSYS' .msh format
- visualization software which can read .vtu file format such as [Paraview](https://www.paraview.org/) or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit), or any tool which can read .vtu file format
- [OpenMPI](https://www.open-mpi.org/) installation (mpif90 for compilation and mpirun for parallel processing)

T-Flows is, in essence, the flow solver without any graphical user interface (GUI).  Although it comes with its own mesh generator, it is very rudimentary and an external software, either free or commercial, would be highly desirable for meshing of complex computational domains.  We regularly use GMSH and would highly recommend it for its inherent scipting ability, but if you have access to any commercial mesh generator which can export meshes in ANSYS' .msh (and .cas, this should be checked) format, that would just fine.  Having no GUI, T-Flows relies on external tools for visualisation of results.  The results are saved in .vtu, Paraview's unstructured data format, and any visualisation software which can read that format is highly desirable for post-processing of results.

From its beginnings, T-Flows was developed for parallel execution with Message Passing Interface (MPI).  If you inted to run it on parallel computational platforms, you will also need an installation of OpenMPI on your system.
<br />

## Optional software packages <a name="soft_req_opt"></a>
<br />

The following packages are not essential to T-Flows, but could prove to be very useful if you become and experienced user, or even developer:

- [grace](https://plasma-gate.weizmann.ac.il/Grace/)
- [PETSc](https://petsc.org/release/)
- git

Visualization tools such as ParaView and VisIt are powerful, self-contained and sufficient for all sorts of post-processings, occasionally you might want to extract profiles from your solution fields and compare them agains experiments or direct numerical simulation (DNS) results, so a two-dimensional plotting software might come handy.  We find grace light particularly suitable for that purpose and many test cases which come with T-Flows, already have benchmark cases compared in xmgrace's format.

Although T-Flows comes with its own suite of linear solvers based on Krylov sub-space family of methods (Incomplete Cholesky and Jacobi preonditioned CG, BiCG and  CGS), to have a better scaling with problem size, you may want to have more choice or even use algebraic multigrid preconditioners available through PETSc.  If PETSc is available on your system, T-Flows' makefiles will link with them and you will have all PETSc solvers at your disposal.

T-Flows resides on [GitHub](www.github.com) platform, and its development is controled by git commands.  Although you can download T-Flows from GitHub as a tarball and use it locally from there on, the connection to GitHub repository gives you the possibility to _pull_ updates, report issues, _track_ your own developments, and even share with them rest of community by pushing your changes.
<br />

# Obtaining the code <a name="obtaining"></a>
<br />

# Compiling the code <a name="compiling"></a>
<br />

