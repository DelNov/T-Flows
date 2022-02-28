# T-Flows 

1. [Introduction](#intro)
2. [Software requirements](#soft_req)
    1. [Minimum](#soft_req_min)
    2. [Higlhy desirable](#soft_req_des)
    3. [Optional](#soft_req_opt)
3. [User requirements](#user_req)
    1. [Minimum](#user_req_min)
    2. [Desirable](#user_req_des)
4. [Obtaining the code](#obtaining)
5. [Compiling the code](#compiling)
    1. [Directory structure](#compiling_dir_struct)
    2. [Sub-programs](#compiling_sub_progs)
6. [Test cases](#test_cases)
    1. [Lid-driven cavity flow](#test_cases_lid_driven_cavity)
    2. [Thermally-driven cavity flow](#test_cases_lid_driven_cavity)
7. [Parallel processing](#parallel_proc)


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

- [git](https://git-scm.com/)
- [PETSc](https://petsc.org/release/)
- [grace](https://plasma-gate.weizmann.ac.il/Grace/)

T-Flows resides on [GitHub](www.github.com) platform, and its development is controled by git commands.  Although you can download T-Flows from GitHub as a tarball and use it locally from there on, the connection to GitHub repository gives you the possibility to _pull_ updates, report issues, _track_ your own developments, and even share with them rest of community by pushing your changes.

Although T-Flows comes with its own suite of linear solvers based on Krylov sub-space family of methods (Incomplete Cholesky and Jacobi preonditioned CG, BiCG and  CGS), to have a better scaling with problem size, you may want to have more choice or even use algebraic multigrid preconditioners available through PETSc.  If PETSc is available on your system, T-Flows' makefiles will link with them and you will have all PETSc solvers at your disposal.

Visualization tools such as ParaView and VisIt are powerful, self-contained and sufficient for all sorts of post-processings, occasionally you might want to extract profiles from your solution fields and compare them agains experiments or direct numerical simulation (DNS) results, so a two-dimensional plotting software might come handy.  We find grace light particularly suitable for that purpose and many test cases which come with T-Flows, already have benchmark cases compared in xmgrace's format.

# User requirements <a name="user_req"></a>

There is no point in denying that successful use of a software package depends on the quality, robustness and intitivity of the software itself, but also on the level of user's experience.  This is even more true for open source software, particularly open source scientific software such as T-Flows.  So, in this section we list some background knowledge required from you in order to successfuly use T-Flows.

## Minimum user requirements <a name="user_req_min"></a>

The bare minimum you need to get T-Flows running is:

- you are able to download T-Flows sources as a .zip file from [GitHub](https://github.com/DelNov/T-Flows)
- you can find your way through your operating system from a terminal, you can find your path to downloaded sources and decompress them
- you know that your operating system has ```make``` command, Fortran and C compiler, and if not, you know who to ask to install them for you

If you are not even on this level, T-Flows is not for you and you are just wasting your time with it.  You are better of venturing into CFD with some commercial package featuring fully fledged GUI.

## Desirable user requirements <a name="user_req_des"></a>

Just like in software requirements section, the minimum will only get you so far.  In order to take a better advantage of T-Flows, your background knowledge should also entail:

- prudence in using Linux operating system from a terminal
- understanding of the ```make``` command
- ability to install third-party software on your computer, such as GMSH and Paraview
- one of the high-level programming languages such as C/C++, Fortran 2003+, Python, Julia or alike
- understanding of fluid mechanics
- essence of finite volume method and how conservation equations are numerically solved and linked
- the basic approaches in turbulence modelling, relative merits of RANS and LES and reasons why hybrid RANS/LES are being used
- approaches to multiphase flows descriptions, in particular VOF for flows with resolved interfaces and Lagrangian particle tracking for flows laden with particles

Ideally, when opting for open-source CFD code such as T-Flows, you should also be:

- eager to see how all components of a CFD program, hence the numerical methods, physical models and linear solvers, are implemented in an unstructured finite volume solver
- able to understand the essence of object-oriented code architecture
- familiar with single-program multiple-data (SPMD) programming paradigm and MPI commands
- version control system with git command
- drive to modify the sources of the program you are using, either in an attempt to improve it or to implement new features be it physical models or numerical methods and algorithms
- ready to adhere to coding standards laid down by T-Flows core development
- willing to share your developments with the rest of the scientific community through GitHub platform

This set of desirable user background knowledge is a bit on the exhaustive side and all the items should not met.  In full openess, not all members of T-Flows core development team are experts in all of these fields, but through organized team work we are covering them all.  Last, but maybe the most import, we hope that the usage of T-Flows will broaden your knowledge in as many of these fields as possible.

# Obtaining the code <a name="obtaining"></a>

The code resides on Github, more preciselly here [T-Flows@Github](https://github.com/DelNov/T-Flows).  You can either download just the zipped sources (from the "Code" dropdown menu, or use git in your terminal to retrive the sources under the git version control system:

    git clone https://github.com/DelNov/T-Flows

> **_Tip 1:_** We strongly recommend the latter approach as git will track all the changes you do to the code and give you access to all possibilities offered by git suite of commands.

If you are using the git command, you might even specify the name of the local directory where you want all the sources to be retrived, for example:

    git clone https://github.com/DelNov/T-Flows  T-Flows-Let-Me-Check-If-It-Works

In any case, the local directory to which all the sources have been retrieved, will be referred to as _root_ from this point in the manual on.

# Compiling the code <a name="compiling"></a>

## Directory structure <a name="compiling_dir_struct"></a>

To cover this section, we assume that you have an open terminal and that you have retreived the sources with one of the two options described in [section above](#obtaining).  The _root_ directory has the following structure:
```
[root]/
├── Binaries
├── Documentation
├── license
├── readme.md
├── Sources
└── Tests
```
> **_Note:_** Remember, [root] is the name of the directory to which you cloned or unzipped the sources.  

The sub-folders have self-explanatory names and we believe that it is only worth mentioning that directory ```Binaries``` will contain executable files.  If you check the contents of the ```Sources``` sub-folder, it will reveal the following structure:
```
[root]/Sources/
├── Convert
├── Divide
├── Generate
├── Libraries
├── Process
├── Shared
└── Utilities
```
which needs a bit more attention.  T-Flows entails four sub-programs called: _Generate_, _Convert_, _Divide_ and _Process_, whose sources lie in corresponding sub-directory names.  Directory ```Shared``` contains sources (by end large Fortran 2008 classes) which are shared by all the programs mentioned above.  Directory ```Libraries``` contains third party libraries and, at the time of writing this manual, contains only [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) libraries for domain decomposition with _Divide_.  Folder ```Utilities``` contains small utilities and prototype procedures used to test some basic concepts and we will epxlain some of them when the need occurs.

## Sub-programs: _Generate_, _Convert_, _Divide_ and _Process_  <a name="compiling_sub_programs"></a>

_Process_, as the name implies, has all the functionality needed to process (discretize and solve) flow equations on a given numerical grid, but it is not a mesh generator on its own.  To provide a mesh to _Process_, you can either use the built-in program _Generate_, which is quite rudimentary and useful for very simple geometries only, or convert a mesh from an external mesh generator (we would recommend [GMSH](https://gmsh.info), but any other software producing the files in Fluent's .msh format will do) using the program _Convert_.  So, as a bare minimum, you have to compile:

1. _Generate_ or _Convert_ and
2. _Process_

Should you want to run your simulations in parallel, you will have to decompose the meshes obtained from _Generate_ or _Convert_ with the program _Divide_.  The parallel processing is not covered here, we dedicate a separate section for it [below](#parallel_proc).  Here, we will only compile the program.

The compilation of each of the sub-programs is performed in their directories (```[root]/Sources/Generate```, ```[root]/Sources/Convert```, ```[root]/Sources/Divide``` and ```[root]/Sources/Process``` by simply issuing command ```make``` in each of them.  When compiling _Convert_, for example, command ```make``` in its directory prints the following on the terminal:
```
#=======================================================================
# Compiling Convert with compiler gnu
#-----------------------------------------------------------------------
# Usage:                                                                
#   make <FORTRAN=gnu/intel/nvidia> <DEBUG=no/yes> OPENMP=<no/yes>      
#        <REAL=double/single>                                           
#                                                                       
# Note: The first item, for each of the options above, is the default.  
#                                                                       
# Examples:                                                             
#   make              - compile with gnu compiler                       
#   make FORTAN=intel - compile with intel compiler                     
#   make DEBUG=yes    - compile with gnu compiler in debug mode         
#   make OPENMP=yes   - compile with gnu compiler for parallel Convert  
#-----------------------------------------------------------------------
gfortran ../Shared/Const_Mod.f90
gfortran ../Shared/Comm_Mod_Seq.f90
gfortran ../Shared/Math_Mod.f90
gfortran ../Shared/File_Mod.f90
...
...
gfortran Load_Fluent.f90
gfortran Load_Gmsh.f90
...
...
gfortran Main_Con.f90
Linking  ../Libraries/Metis_5.1.0_Linux_64/libmetis_i32_r64.a ../../Binaries/Convert 
```
At each invokation of the command ```make``` for any of the four programs` folders, makefile will print a header with possible options which can be passed to ```make```.  In the above you can see that you can specify which compiler you want to use (we use almost exclusivelly gnu), wheather you want debugging information in the executable and, probably most important of all, the precision of the floating point numbers.  If you specify ```make REAL=single```, ```make``` will compile the program with 32-bit representation of floating point numbers and if you define ```make REAL=double```, it will use the 64-bit representation.  

The command ```make clean``` will clean all object and module files from the local directory.

> **_Note 1:_** You don't really have to specify ```REAL=double``` or ```FORTRAN=gnu``` since these are the default options, as described in the header written by the makefile above.

> **_Note 2:_** Although the provided ```makefiles``` take care of the dependency of the sources, if you change precision in-between two compilations, you will have to run a ```make clean``` in between to make sure that objects with 32-bit and 64-bit representation of floating points mix up.

> **_Warning:_** All sub-programs should be compiled with the same precision.  At the time of writing these pages, _Process_ compiled with double precision will not be able to read files created by _Convert_ in single precision.  So, decide on a precision and stick with it for all sub-programs.

After you visit all four sub-directories with sources and run the ```make``` command in each, the sub-folder ```[root]/Binaries``` will have the following files:
```
[root]/Binaries/
├── Convert
├── Divide
├── Generate
├── Process
└── readme
```
and you are ready to start your first simulations.

# Test cases <a name="test_cases"></a>

## Lid-driven cavity flow <a name="test_cases_lid_driven_cavity"></a>

It is hard to imagine a problem in CFD simpler than a lid-driven flow in a cavity:

![Lid-driven cavity domain!](Documentation/Manual/Figures/lid_driven_cavity_domain_050.png "Lid driven cavity domain")

The flow occurs in a cavity with square cross section, in which all the walls are static except the top one which is moving.  The cavity is long enough in the spanwise direction that an assumption of two-dimensionality can be made.  The simplicity of this case is hard to beat as it occurs in one of the simplest possible domain geometries (a cube or a square in two dimensions), the fact that there are no inflow and outflow boundaries, and that there are steady laminar solutions for a range of Reynolds numbers.  Owing to its simplicity, it is well benchmarked and quite often the first test case one ever considers with a new CFD code.

We will use this case to introduce a few new concepts in T-Flows:
- conversion of a grid from GMSH with _Convert_
- specification of periodic boundary conditions
- setting up a control file for _Process_
- running the simulation and adjusting the control file
- visualization of results

In order to run this case, please go to the directory ```[root]/Tests/Manual/Lid_Driven_Cavity```.  There you will find the following files:
* ```lid_driven_cavity.geo```
* ```lid_driven_cavity.msh```

The first file, with extension ```.geo``` is the script GMSH uses to generate a mesh.  It is beyond the scope of this manual to teach you how to use third party software, but it's probably worth having a look at it:
```
A  =  1.0;  // length and height of the cavity
B  =  0.5;  // width of the cavity
NA = 20;    // resolution in length and height
NB =  3;    // resolution in width (periodic direction)

// Define points
Point(1) = {0, 0, 0};  Point(2) = {A, 0, 0};
Point(3) = {A, 0, A};  Point(4) = {0, 0, A};

// Define lines
Line(1) = {1, 2};  Line(2) = {2, 3};
Line(3) = {3, 4};  Line(4) = {4, 1};

// Define front surface
Curve Loop(1) = {1, 2, 3, 4};  Plane Surface(1) = {1};

// Set lines to transfinite
// (+1 is because GMSH expects number of points here)
Transfinite Curve {1, 2, 3, 4} = NA+1  Using Bump  1.0;

// Set surface to be transfinite and quadrilateral
Transfinite Surface {1} = {1, 2, 3, 4};  Recombine Surface {1};

// Expand in spanwise direction
Extrude {0, B, 0} {Surface{1}; Layers {NB}; Recombine;}

// Define boundary conditions
Physical Surface("moving_wall", 27) = {21};
Physical Surface("static_wall", 28) = {25, 13, 17};
Physical Surface("periodic_y", 29) = {1, 26};

// Define volume condition
Physical Volume("fluid", 30) = {1};
```
to see that it is very readable and intuitive.  Script defines the points, the lines, a surface, extends this surface to create a volume in third dimension and finally prescribes the boundary conditions (called ```Pysical Surface``` and set to ```moving_wall```, ```static_wall``` and ```periodicity_y``` for this case.  The script is wrapped up by specifying a ```Physical Volume```, whose name is simply set to ```fluid```.



## Thermally-driven cavity flow <a name="test_cases_thermally_driven_cavity"></a>

# Parallel processing  <a name="parallel_proc"></a>

