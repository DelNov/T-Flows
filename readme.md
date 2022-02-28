# T-Flows 

1. [Introduction](#intro)
2. [Software requirements](#soft_req)
    1. [Minimim](#soft_req_min)
    2. [Higlhy desirable](#soft_req_des)
    3. [Optional](#soft_req_opt)
2. [User requirements](#user_req)
    1. [Minimim](#user_req_min)
    2. [Higlhy desirable](#user_req_des)
3. [Obtaining the code](#obtaining)
4. [Compiling the code](#compiling)


# Introduction <a name="intro"></a>

T-Flows is a computational fluid dynamics (CFD) program for simulation of turbulent, single and multiphase flows.  Numerical method is based on collocated finite volume method on unstructured arbitrary grids and turbulence models include a range of Reynolds-averaged Navier-Stokes (RANS) models, large eddy simulations (LES), as well as hybrid RANS-LES approach.  A more comprehensive list of turbulence models is [here](https://github.com/DelNov/T-Flows/blob/bojan_petsc_solvers_almost_alpha/Documentation/Manual/turbulence_models.md).

Multiphase models include an algebraic volume of fluid (VOF) method and Lagrangian particle tracking model.  Three-phase flows situations (two fluid phases with VOF and one solid phase as particles) are also supported.

> **_Note:_** In T-Flows, the Navier-Stokes equations are discretized in their _incompressible_ form, meaning _only_ that pressure and temperatures are _not_ linked through an equation of state.  All physical properties in T-Flows can be variable, but you should keep in mind that variable density does not mean compressibility.

# Software requirements <a name="soft_req"></a>

## Minimum software requirements <a name="soft_req_min"></a>

The bare minimum to get T-Flows running entails:

- make utility
- Fortran 2008 compiler
- standard C compiler

T-Flows is almost entirely written in Fortran 2008 (only one function is written in C) and the compilation is controlled by makefiles.  So, the the requirements listed above are a bare minimum for you to start using the code.  
t
Although there is, in principle, no restriction on the operating system on which you can use T-Flows, its natural habitat is Linux, as we develop test it on Linux, and Linux meets the minimum software requirements either _out of the box_, or with minimum installation effort.

> **_Note:_** We do not specify the minimum version for any of the required or recommended software.  We believe that if you are reading these pages, you do have access to a relatively recent hardware which also implies an up to date operating system and the associated tools.  

## Highly desirable software requirements <a name="soft_req_des"></a>

Although without meeting the minimum software requirements listed above you will not get anywhere, they alone will not get you very far either.  To make a practical use of T-Flows, it is highly desirable that you also have the following:

- [GMSH](https://gmsh.info)
- any other free or commercial mesh generator exporting ANSYS' .msh format
- visualization software which can read .vtu file format such as [Paraview](https://www.paraview.org/) or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit), or any tool which can read .vtu file format
- [OpenMPI](https://www.open-mpi.org/) installation (mpif90 for compilation and mpirun for parallel processing)

T-Flows is, in essence, the flow solver without any graphical user interface (GUI).  Although it comes with its own mesh generator, it is very rudimentary and an external software, either free or commercial, would be highly desirable for meshing of complex computational domains.  We regularly use GMSH and would highly recommend it for its inherent scipting ability, but if you have access to any commercial mesh generator which can export meshes in ANSYS' .msh (and .cas, this should be checked) format, that would just fine.  Having no GUI, T-Flows relies on external tools for visualisation of results.  The results are saved in .vtu, Paraview's unstructured data format, and any visualisation software which can read that format is highly desirable for post-processing of results.

From its beginnings, T-Flows was developed for parallel execution with Message Passing Interface (MPI).  If you inted to run it on parallel computational platforms, you will also need an installation of OpenMPI on your system.


## Optional software packages <a name="soft_req_opt"></a>

The following packages are not essential to T-Flows, but could prove to be very useful if you become and experienced user, or even developer:

- [grace](https://plasma-gate.weizmann.ac.il/Grace/)
- [PETSc](https://petsc.org/release/)
- git

Visualization tools such as ParaView and VisIt are powerful, self-contained and sufficient for all sorts of post-processings, occasionally you might want to extract profiles from your solution fields and compare them agains experiments or direct numerical simulation (DNS) results, so a two-dimensional plotting software might come handy.  We find grace light particularly suitable for that purpose and many test cases which come with T-Flows, already have benchmark cases compared in xmgrace's format.

Although T-Flows comes with its own suite of linear solvers based on Krylov sub-space family of methods (Incomplete Cholesky and Jacobi preonditioned CG, BiCG and  CGS), to have a better scaling with problem size, you may want to have more choice or even use algebraic multigrid preconditioners available through PETSc.  If PETSc is available on your system, T-Flows' makefiles will link with them and you will have all PETSc solvers at your disposal.

T-Flows resides on [GitHub](www.github.com) platform, and its development is controled by git commands.  Although you can download T-Flows from GitHub as a tarball and use it locally from there on, the connection to GitHub repository gives you the possibility to _pull_ updates, report issues, _track_ your own developments, and even share with them rest of community by pushing your changes.


# User requirements <a name="user_req"></a>

There is no point in denying that successful use of a software package depends on the quality, robustness and intitivity of the software itself, but also on the level of user's experience.  This is even more true for open source software, particularly open source scientific software such as T-Flows.  So, in this section we list some background knowledge required from you in order to successfuly use T-Flows.

## Minimum user requirements <a name="user_req_min"></a>

The bare minimum you need to get T-Flows running is:

- you can use a web browser and download a tarball from GitHub
- you can find your way through your operating system from a terminal, you are able to create directories, change them, find downloaded files and unpack them if they are compressed
- you know that your operating system has make command, Fortran and C compiler, and if not, you know who to ask to install them for you.

If you are not on this level, T-Flows is not for you and you are just wasting your time with it.

## Highly desirable user requirements <a name="user_req_des"></a>

Just like in software requirements section, the minimum will only get you so far.  In order to take a better advantage of T-Flows, your background knowledge also entails:

- the operating system on which you prudently manouver happens to be a Linux
- you understand that using a keyboard and a terminal is much faster than using a mouse and windows
- you understand how make command works, how is the makefile organized, and that you can pass options to makefile
- you are able to install missing software on the computer you want to run T-Flows on
- you have a basic understanding of one of the high-level programming languages such as C/C++, Fortran 2003+, Python, Julia or alike
- you have basic understanding of fluid mechanics
- you know the essence of finite volume method and how conservation equations are numerically solved and linked
- you know that turbulence modelling is an important part of CFD and you've heard about approaches such as RANS and LES
- you know that there are many approaches to simulate multiphase flows and that VOF simulates interfaces between phases, whereas Lagrangian particle tracking is used to track individual particles in the carrier fluid.

## Perfect user requirements <a name="user_req_per"></a>

You will take the best advantage of T-Flows if, in addition to the requirements listed so far, the following is also true:

- you are eager to see how all components of a CFD program, hence the numerical methods, physical models and linear solvers, are implemented in an unstructured finite volume solver
- you are not a stranger to an object-oriented architecture of the code
- you are familiar with single-program multiple-data (SPMD) programming paradigm and understand MPI commands
- you are eager to modify the program you are using, either in an attempt to improve it or implement new features be it physical models or numerical algorithms
- you underestand that deviationg from coding standards laid down by T-Flows core development team is _not_ creativity.
- you are ready to share your developments with the rest of the scientific community through GitHub platform
- you believe that counting from one is more intuitive than counting from zero
- you are aware that Fortran is passing arguments by reference, but you think this notion is better than constant contemplations at each function call if you should pass by value, by pointer or by reference.
- 

# Obtaining the code <a name="obtaining"></a>


# Compiling the code <a name="compiling"></a>


