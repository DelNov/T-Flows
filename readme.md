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
    2. [Sub-programs](#compiling_sub_programs)
6. [Demonstration cases](#test_cases)
    1. [Lid-driven cavity flow](#demo_lid_driven)
        1. [On hexahedral grid](#demo_lid_driven_hexa)
            1. [Creating the grid](#demo_lid_driven_hexa_create)
            1. [Converting the grid](#demo_lid_driven_hexa_convert)
        2. [On polyhedral grid](#demo_lid_driven_dual)
    2. [Thermally-driven cavity flow](#demo_thermally_driven)
        1. [With Boussinesq approximation](#demo_thermally_driven_boussinesq)
        2. [With variable physical properties](#demo_thermally_driven_variable)
    3. [Parallel processing](#demo_parallel_proc)
        1. [Compiling for parallel runs](#demo_parallel_proc_compiling)
        2. [Creating and dividing the grid](#demo_parallel_proc_dividing)
        3. [Running the simulation in parallel](#demo_parallel_proc_running)
7. [Benchmark cases](#bench_cases)
    1. [Large eddy simulation over a matrix of cubes](#bench_cases_matrix)


# Introduction <a name="intro"></a>

T-Flows is a computational fluid dynamics (CFD) program for simulation of turbulent,
single and multiphase flows.  Numerical method is based on collocated finite
volume method on unstructured arbitrary grids and turbulence models include a
range of Reynolds-averaged Navier-Stokes (RANS) models, large eddy simulations
(LES), as well as hybrid RANS-LES approach.  A more comprehensive list of
turbulence models is [here](https://github.com/DelNov/T-Flows/blob/bojan_petsc_solvers_almost_alpha/Documentation/Manual/turbulence_models.md).

Multiphase models include an algebraic volume of fluid (VOF) method and
Lagrangian particle tracking model.  Three-phase flows situations (two fluid
phases with VOF and one solid phase as particles) are also supported.

> **_Note:_** In T-Flows, the Navier-Stokes equations are discretized in their
_incompressible_ form, meaning _only_ that pressure and temperatures are _not_
linked through an equation of state.  All physical properties in T-Flows can
be variable, but you should keep in mind that variable density does not mean
compressibility.

# Software requirements <a name="soft_req"></a>

## Minimum software requirements <a name="soft_req_min"></a>

The bare minimum to get T-Flows running entails:

- make utility
- Fortran 2008 compiler
- standard C compiler

T-Flows is almost entirely written in Fortran 2008 (only one function is written
in C) and the compilation is controlled by makefiles.  So, the the requirements
listed above are a bare minimum for you to start using the code.  

Although there is, in principle, no restriction on the operating system on
which you can use T-Flows, its natural habitat is Linux.   We develop test T-Flows
on Linux, and Linux meets the minimum software requirements either _out of the box_,
or with small installation effort.

> **_Note:_** We do not specify the minimum version for any of the required or
recommended software.  We believe that if you are reading these pages, you do
have access to a relatively recent hardware which also implies an up to date
operating system and the associated tools.  

## Highly desirable software requirements <a name="soft_req_des"></a>

Although without meeting the minimum software requirements listed above you
will not get anywhere, they alone will not get you very far either.  To make a
practical use of T-Flows, it is highly desirable that you also have the following:

- [GMSH](https://gmsh.info)
- any other free or commercial mesh generator exporting ANSYS' ```.msh``` or
  ```.cas``` format (**double-check** for .cas)
- visualization software which can read ```.vtu``` file format such as
  [ParaView](https://www.paraview.org/) or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit),
  or any tool which can read ```.vtu``` file format
- [OpenMPI](https://www.open-mpi.org/) installation (mpif90 for compilation and
  mpirun for parallel processing)

T-Flows is, in essence, the flow solver without any graphical user interface (GUI).  
Although it comes with its own mesh generator, it is very rudimentary and an
external grid generation software would be highly desirable for meshing complex
computational domains.  We regularly use GMSH and would highly recommend it for
its inherent scipting ability, but if you have access to any commercial grid
generator which can export meshes in ANSYS' ```.msh``` (and .cas, this should
be checked) format, that would just fine.  Having no GUI, T-Flows relies on
external tools for visualisation of results.  The results are saved in ```.vtu```,
ParaView's unstructured data format, and any visualisation software which can
read that format is highly desirable for post-processing of results.

From its beginnings, T-Flows was developed for parallel execution with Message
Passing Interface (MPI).  If you inted to run it on parallel computational
platforms, you will also need an installation of OpenMPI on your system.

## Optional software packages <a name="soft_req_opt"></a>

The following packages are not essential to T-Flows, but could prove to be very
useful if you become and experienced user, or even developer:

- [git](https://git-scm.com/)
- [PETSc](https://petsc.org/release/)
- [grace](https://plasma-gate.weizmann.ac.il/Grace/)

T-Flows resides on [GitHub](www.github.com) platform, and its development is
controlled by ```git``` commands.  Although you can download T-Flows from
GitHub as a ```.zip``` file and use it locally from there on, the connection
to GitHub repository gives you the possibility to _pull_ updates, report issues,
_track_ your own developments, and even share with them rest of community by
pushing your changes.

Although T-Flows comes with its own suite of linear solvers based on Krylov
sub-space family of methods (Incomplete Cholesky and Jacobi preconditioned
conjugate gradient (CG), bi-conjugate gradient (BiCG) and conjugate gradeint
squared (CGS)), you might want to have more choice or even want to use
algebraic multigrid preconditioners for better scaling with problem size.  
To that end, we linked T-Flows with PETSc.  If PETSc is installed on your
system, T-Flows' makefiles will recognize, link with PETSc libraries and you
will have all PETSc solvers at your disposal.

Visualization tools such as ParaView and VisIt are powerful, self-contained
and sufficient for all sorts of post-processing, occasionally you might want
to extract profiles from your solution fields and compare them against
experiments or direct numerical simulation (DNS) results, so a two-dimensional
plotting software might come handy.  We find grace light particularly suitable
for that purpose and many test cases which come with T-Flows, already have
benchmark cases compared in grace's format (```.agr``` files).

> **_Note:_** During its lifespan, grace went through a number of transitions
and is known under different names.  It used to be called _xmgr_, _xmgrace_,
_ace_ and currently _grace_.  It is possible that, with installation of grace,
its other names become available on your system as aliases.

# User requirements <a name="user_req"></a>

There is no point in denying that successful use of a software package depends
on the quality, robustness and intuitivity of the software itself, but also on
the user's previous experience.  This is even more true for open source
software, particularly open source scientific software such as T-Flows.  
So, in this section we list some background knowledge required from you in
order to successfully use T-Flows.

## Minimum user requirements <a name="user_req_min"></a>

The bare minimum you need to get T-Flows running is:

- you are able to download T-Flows sources as a .zip file from [GitHub](https://github.com/DelNov/T-Flows)
- you don't need to click icons in a file manager to reach your files or execute
  commands, because you can use your operating system (OS) from a terminal.  
- you know that your operating system has ```make``` command, Fortran and
  C compiler, and if not, you know who to ask to install these facilites for you

If you are not even on this level, T-Flows is not for you and you are just
wasting your time with it.  You are better of venturing into CFD with some
commercial package featuring fully fledged GUI.

## Desirable user requirements <a name="user_req_des"></a>

Just like in software requirements section, the minimum will only get you so far.  
In order to take a better advantage of T-Flows, your background knowledge should
also entail as many of the following:

- prudence in using Linux operating system from a terminal
- understanding of the ```make``` command
- ability to install third-party software on your computer, such as GMSH, ParaView
  or maybe even build PETSc from the sources
- one of the high-level programming languages such as C/C++, Fortran 2003+,
  Python, Julia or alike
- understanding of fluid mechanics
- essence of finite volume method and how conservation equations are numerically
  solved and linked
- the basic approaches in turbulence modeling, relative merits of RANS and LES
  and reasons why hybrid RANS/LES are being used
- approaches to multiphase flows description, in particular VOF for flows with
  resolved interfaces and Lagrangian particle tracking for flows laden with particles

Ideally, when opting for open-source CFD code such as T-Flows, you should also be:

- eager to see how all components of a CFD program, hence the numerical methods,
  physical models and linear solvers, are implemented in an unstructured finite
  volume solver
- able to understand the essence of object-oriented code architecture
- familiar with single-program multiple-data (SPMD) programming paradigm and
  related MPI commands
- version control system with ```git``` command
- drive to modify the sources of the program you are using, either in an
  attempt to improve it or to implement new features be it physical models or
  numerical methods and algorithms
- ready to adhere to coding standards laid down by T-Flows core development team
- willing to share your developments with the rest of the scientific community
  through GitHub platform

This set of desirable user background knowledge is a bit on the exhaustive side
and all the items should not met.  In full openness, not even all members of
T-Flows core development team are experts in _all_ of these fields, but through
organized team work we are covering them all.  Last, but maybe the most import,
we hope that the usage of T-Flows will broaden your knowledge in as many of these
fields as possible, and that you will find this journey exciting.

# Obtaining the code <a name="obtaining"></a>

The code resides on Github, more precisely here: [T-Flows](https://github.com/DelNov/T-Flows).  
You can either download just the zipped sources (from the "Code" drop-down menu,
or use git in your terminal to retrieve the sources under the git version control
system:
```
git clone https://github.com/DelNov/T-Flows
```

> **_Tip 1:_** We strongly recommend the latter approach as ```git``` will track all
the changes you do to the code and give you access to all possibilities offered
by git suite of commands.

If you are using the ```git``` command, you might even specify the name of the
local directory where you want all the sources to be retrieved, for example:
```
git clone https://github.com/DelNov/T-Flows  T-Flows-Let-Me-Check-If-It-Works
```

In any case, the local directory to which all the sources have been retrieved,
will be referred to as ```[root]``` from this point in the manual on.

# Compiling the code <a name="compiling"></a>

## Directory structure <a name="compiling_dir_struct"></a>

To cover this section, we assume that you have an open terminal and that you
have retrieved the sources with one of the two options described in
[section above](#obtaining).  The ```[root]``` directory has the following
structure:
```
[root]/
├── Binaries
├── Documentation
├── license
├── readme.md
├── Sources
└── Tests
```

> **_Note:_** Remember, [root] is the name of the directory to which you cloned
or unzipped the sources.  It is not ```/``` in Linux or ```C:\``` on Windows.

The sub-folders have reather self-explanatory names and we believe that it is
only worth mentioning that the directory ```Binaries``` will contain executable
files.  If you check the contents of the ```Sources``` sub-folder, it will reveal
the following structure:
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

which needs a bit more attention.  T-Flows entails four sub-programs called:
_Generate_, _Convert_, _Divide_ and _Process_, whose sources lie in corresponding
sub-folder names.  Folder ```Shared``` contains sources (by end large Fortran
2008 _classes_) which are shared by all the programs mentioned above.  
Folder ```Libraries``` contains third party libraries and, at the time of writing
this manual, contains only [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
libraries for domain decomposition with _Divide_.  Folder ```Utilities``` contains
small utilities and prototype procedures used to test some basic concepts or help
with some menial takss and we will explain some of them when the need occurs.

## Sub-programs: _Generate_, _Convert_, _Divide_ and _Process_  <a name="compiling_sub_programs"></a>

_Process_, as the name implies, has all the functionality needed to process
(discretize and solve) flow equations on a given numerical grid, but it is not
a grid generator on its own.  To provide a grid to _Process_, you can either
use the built-in program _Generate_, which is quite rudimentary and useful for
very simple geometries only, or convert a mesh from an external mesh generator.

> **_Note:_**  As third-party grid generator we would recommend [GMSH](https://gmsh.info),
although any other software producing the files in Fluent's ```.msh``` format will do)

So, as a bare minimum, you have to compile:

1. _Generate_ or _Convert_ and
2. _Process_

Should you want to run your simulations in parallel, you will have to decompose
the grids obtained from _Generate_ or _Convert_ with the program _Divide_.  
The parallel processing is not covered here, we dedicate a separate section
for it [below](#demo_parallel_proc).  Here, we will only compile the program _Divide_.

The compilation of each of the sub-programs is performed in their directories
(```[root]/Sources/Generate```, ```[root]/Sources/Convert```, ```[root]/Sources/Divide```
and ```[root]/Sources/Process``` by simply issuing command ```make``` in each
of them.  When compiling _Convert_, for example, command ```make``` in its
directory prints the following on the terminal:
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
At each invokation of the command ```make``` for any of the four program's
folders, ```makefile``` will print a header with possible options which can be
passed to ```make```.  In the above you can see that you can specify which
compiler you want to use (we use almost exclusively gnu), whether you want
debugging information in the executable and, probably most important of all,
the precision of the floating point numbers.  If you specify ```make REAL=single```,
```make``` will compile the program with 32-bit representation of floating
point numbers and if you define ```make REAL=double```, it will use the 64-bit
representation.  

The command ```make clean``` will clean all object and module files from the local directory.

> **_Note 1:_** You don't really have to specify ```REAL=double``` or
```FORTRAN=gnu``` since these are the default options, as described in the header
written by the ```makefile``` above.

> **_Note 2:_** Although the provided ```makefiles``` take care of the dependencies
of the sources, if you change precision in-between two compilations, you will have
to run a ```make clean``` in between to make sure that objects with 32-bit and
64-bit representation of floating points mix up.  

> **_Warning:_** All sub-programs should be compiled with the same precision.  
At the time of writing these pages, _Process_ compiled with double precision
will not be able to read files created by _Convert_ in single precision.  So,
decide on a precision and stick with it for all sub-programs.

After you visit all four sub-directories with sources and run the ```make```
command in each, the sub-folder ```[root]/Binaries``` will contain the following
files:
```
[root]/Binaries/
├── Convert
├── Divide
├── Generate
├── Process
└── readme
```
and you are ready to start your first simulations.

# Demonstration cases <a name="test_cases"></a>

The cases presented in this chapter serve only to show the basic usage of T-Flows.
They are by no means adhering to best practice guidlines, nor is their aim to be
used as verifcation against experiments, DNS or other published results.  
The following chapter is dedicated to benchmarking and more rigirous CFD
practices.

## Lid-driven cavity flow <a name="demo_lid_driven"></a>

It is hard to imagine a problem in CFD simpler than a lid-driven flow in a
cavity:

![Lid-driven cavity domain!](Documentation/Manual/Figures/lid_driven_domain.png "Lid-driven cavity domain")

The flow occurs in a cavity with square cross section, in which all the walls
are static except the top one which is moving.  The cavity is long enough in
the spanwise direction that an assumption of two-dimensionality can be made.  
The simplicity of this case is hard to beat as it occurs in one of the simplest
possible domain geometries (a cube or a square in two dimensions), the fact
that there are no inflow and outflow boundaries, and that there are steady
laminar solutions for a range of Reynolds numbers.  Owing to its simplicity,
it is well documented, often used to benchmark CFD codes and quite often the
first test case one ever considers with a new CFD code.

We will use this case to introduce a few new concepts in T-Flows:
- conversion of a grid from GMSH with _Convert_
- specification of periodic boundary conditions
- setting up a control file for _Process_
- running the simulation and adjusting the control file
- visualization of results
- typical workflow in T-Flows

### On hexahedral grid <a name="demo_lid_driven_hexa"></a>

#### Creating the mesh <a name="demo_lid_driven_hexa_create"></a>

We start with a grid build from hexahedral cell. In order to run this case,
please go to the directory ```[root]/Tests/Manual/Lid_Driven_Cavity/Hexa```.  
There you will find the following files:
```
Hexa/
├── lid_driven.geo
└── lid_driven_msh.gz
```

The first file, with the extension ```.geo```, is the script GMSH uses to
generate a mesh, whereas the second file, ```lid_driven.msh.gz```, is in this
directory in case you don't have GMSH installed on your system.  If that is the
case, skip directly [here](#demo_lid_driven_hexa_convert).

Although it is beyond the scope of this manual to teach you how to use third
party software, we believe it's interesting to have a look at the
```lid_driven.geo``` file:
```
A  =  1.0;  // length and height of the cavity
B  =  0.1;  // width of the cavity
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
Physical Surface("upper_wall", 27) = {21};
Physical Surface("side_walls", 28) = {25, 17};
Physical Surface("lower_wall", 29) = {13};
Physical Surface("periodic_y", 30) = {1, 26};

// Define volume condition
Physical Volume("fluid", 31) = {1};

```

to see that it is very readable and intuitive.  Script defines the points,
the lines, a surface, extends this surface to create a volume in third
dimension and finally prescribes the boundary conditions (called ```Pysical Surface```
and set to ```upper_wall```, ```side_walls```, ```lower_wall``` and ```periodicity_y```
for this case.  The script is wrapped up by specifying a ```Physical Volume```,
whose name is simply set to ```fluid```.

> **_Note:_** GMSH comes with a GUI, and you don't have to write a script like
that at all - you can create a mesh solely from GUI.  That might seem comfortable
at first, but it comes with a big disadvantage, particularly for scientific work,
that it is _not_ reproducible.  However, the main strength of GMSH as we see it
is that it automatically creates a script with commands you enter in GUI so that
you can use it for reproducing your actions, or tuning the script for more
precise adjustments.  We usually create grids in GMSH in an iterative fashion:
work in the GUI for a few steps, then exit and examine the script it created or
updated and make fine tuning as we please.  Then we switch back to GUI with
updated script and make more steps there.  Another advantage of GMSH's scripting
is the possibility of parametrization.  We introduced variables ```A```, ```B```,
```NA``` and ```NB``` to be able to easily change domain dimensions and resolutions.

> **_Warning:_** It is easy to forget to specify the ```Physical Volume```
in GMSH which is a pitfall, because if you don't do it, it will only export
surface mesh blocking T-Flows in discretizing equations through volumes.

Since the script for GMSH is ready, you can run it from command line with this command:
```
gmsh -3 lid_driven.geo
```
which will, after a few instances, create the file:
```
lid_driven.msh
```
> **_Note:_** Both GMSH and Fluent produce grid files with extnesion ```.msh```,
the formats are completely different and GMSH makes an educated guess which one
it is reading based on its contents.

#### Converting the mesh to T-Flows format <a name="demo_lid_driven_hexa_convert"></a>

Once you have the mesh in this format, you can use _Convert_ to convert it
into T-Flows native file format.  This will take some and will require your
attention, so we don't recommend you going on with _Convert_ unless you are
sure you can have 15 minutes of peace and focus ahead of you.  We assume you
compiled _Convert_ as it was described in [this section](#compiling_sub_programs).   
Given that you are in directory ```[root]/Tests/Manual/Lid_Driven_Cavity/Hexa```,
you can call _Convert_ with:
```
../../../../Binaries/Convert
```

which will prompt you with the screen:
```
 #=======================================================================
 #                                                                   
 #    ______________________ ____    ________  __      __  _________
 #    \__    ___/\_   _____/|    |   \_____  \/  \    /  \/   _____/
 #      |    |    |    __)  |    |    /   |   \   \/\/   /\_____  \  
 #      |    |    |     \   |    |___/    |    \        / /        \
 #      |____|    \___  /   |_______ \_______  /\__/\  / /_______  /
 #                    \/            \/       \/      \/          \/  
 #                     _____                      __                 
 #                    / ___/__  ___ _  _____ ____/ /_                
 #                   / /__/ _ \/ _ \ |/ / -_) __/ __/                
 #                   \___/\___/_//_/___/\__/_/  \__/                 
 #                                                                   
 #                         Double precision mode
 #-----------------------------------------------------------------------
 #================================================================
 # Enter the name of the grid file you are importing (with ext.):
 #----------------------------------------------------------------
```

Here you have to type the name of the grid with extension, hence ```lid_driven.msh```.  
Many messages on the screen will follow outlining the conversion process which we
won't cover now, but the _Convert_ will stop with the next question:
```
 #=================================================
 # Would you like to create a dual grid? (yes/no)
 #-------------------------------------------------
```

at which point you type ```no```.  Dual grids will be covered in the
[next example](#demo_lid_driven_dual).  Next question _Convert_ asks you
concerns geometric extents:
```
 #=========================================
 # Geometric extents:                 
 #-----------------------------------------
 # X from:  0.000E+00  to:  1.000E+00
 # Y from:  0.000E+00  to:  5.000E-01
 # Z from:  0.000E+00  to:  1.000E+00
 # Enter scaling factor for geometry:
 # (or skip to keep as is):
 #-----------------------------------------
```

You may wish to scale the entire geometry if, for example, you generated it in
millimeters but you know that T-Flows' _Process_ works in SI units, or if you
want to use scaling to change the non-dimensional numbers (such as Reynolds or
Rayleigh) that way.  For this case, we are not doing it, so just answer ```skip```.  
Next question which _Convert_ will pose you concerns connections of boundary
cells to inside cells:
```
 #====================================
 # Position the boundary cell centres:
 #------------------------------------
 # Type 1 for barycentric placement
 # Type 2 for orthogonal placement
 #------------------------------------
```

We noticed that, in some cases, the accuracy of gradient computations with the
least squares method is more accurate if these connections are orthogonal, so
it is safe to answer with ```2```.

> **_Note:_** For orthogonal grids, these boundary cell centers coincide, so it
hardly matters.  It makes a difference on distorted grids.

Next thing _Convert_ will ask you is very important and is about the periodic
boundary conditions:
```
 #======================================================
 # Grid currently has the following boundary conditions:
 #------------------------------------------------------
 #  1. UPPER_WALL
 #  2. SIDE_WALLS
 #  3. LOWER_WALL
 #  4. PERIODIC_Y
 #------------------------------------------------------
 #==============================================================
 # Enter the ordinal number(s) of periodic-boundary condition(s)
 # from the boundary condition list (see above)                 
 # Type skip if there is none !                                 
 #--------------------------------------------------------------
```

It lists all the boundary conditions specified in the ```.msh``` file and asks
you which of those you want to set to be periodic.  In this case, you enter ```4```,
because the periodicity set in GMSH to be called ```periodic_y``` is fourth in
the list.  Since you can have periodicity in more than one direction, _Convert_
will prompt you with the same question but with a shorter list of boundary conditions:
```
 #======================================================
 # Grid currently has the following boundary conditions:
 #------------------------------------------------------
 #  1. UPPER_WALL
 #  2. SIDE_WALLS
 #  3. LOWER_WALL
 #------------------------------------------------------
 #==============================================================
 # Enter the ordinal number(s) of periodic-boundary condition(s)
 # from the boundary condition list (see above)                 
 # Type skip if there is none !                                 
 #--------------------------------------------------------------
```
and this time you answer ```skip```.  Finally, _Convert_ asks you about the
distance from the wall calculation:

```
 #==================================================================
 # Calculating distance from the walls                              
 #------------------------------------------------------------------
 # Type ordinal number(s) of wall or wall_flux boundary condition(s)
 # from the boundary condition list (see above) separated by space.
 # Cells' centers distances to the nearest wall will be calculated
 # for the listed wall boundary(s).                                 
 #                                                                  
 # This is needed for RANS and HYBRID turbulence models as well as  
 # for proper initialization with potential pressure-like field.    
 #                                                                  
 # Type skip to skip this and set wall distance to one everywhere.  
 #------------------------------------------------------------------
```
Distance to the wall is important for many turbulence models as well as for
problems with conjugate heat transfer.  It can be computed in two ways:
- Geometrically from _Convert_, browsing through all inside cells and searching
  for nearest boundary cells which is accurate but can be time consuming for
  large grids, or
- From a partial differential equation proposed by [Spalding](https://www.cfd-online.com/Wiki/Transport_equation_based_wall_distance_calculation)
  which is done from _Process_ in cases the above step is skipped in _Convert_.  
  Wall distance calculation from partial differential equation is faster for
  large grids, but less accurate.  However, the errors in the near wall regions
  are rather small, invariantly smaller than 1%.

For this case, let's compute them from here so you answer ```1 2``` to the
above question to instruct _Convert_ which boundaries can be considered as
solid walls.  

Before it ends, _Convert_ will ask you one more thing:
```
 #===========================================
 # Creating 1d file with the node
 # coordinates in non-homogeneous directions
 #-------------------------------------------
 # Insert non-homogeneous direction
 # (x, y, z, rx, ry, rz or skip)
 # -------------------------------------------
```

which is important for computation of turbulent flows and described with some
turbulent benchmark cases below.  For the time being, feel free to answer with
a ```skip```.

#### Analyzing the results of the _Convert_

During the conversion process, _Convert_ creates the following files:
```
Hexa/
├── lid_driven.cfn
├── lid_driven.dim
├── lid_driven.vtu
├── lid_driven.faces.vtu
├── lid_driven.shadows.vtu
└── control_template_for_lid_driven
```

Which will be described next.  Files with extension ```.cfn``` and ```.dim```
are binary files in internal T-Flows format, which can be read by _Process_
(or _Divide_ for domain decomposition, which is covered in the separate
[section](#demo_parallel_proc) below).  A number of files with extension ```.vtu```
is created too.  You can visualize them with ParaView or VisIt for inspection.  
The most interesting is the ```lid_driven.faces.vtu``` because it shows
boundary conditions.  Once visualized, the ```lid_driven.faces.vtu``` shows
the following:

![Lid-driven hexa front!](Documentation/Manual/Figures/lid_driven_hexa_front.png "Lid-driven hexa front")

which reveals boundary conditions for this case, shown here in different colors.  
Light blue is for the inside cells (no boundary conditions there) and orange
and red are for different boundary conditions.

if you rotate the domain in ParaView, you will see something which may surprise
you at first:

![Lid-driven hexa back!](Documentation/Manual/Figures/lid_driven_hexa_back.png "Lid-driven hexa back")

The cells in the back seem to be hollow, as if they are missing the faces on
periodic boundary.  This is done on purpose.  Since T-Flows uses face-base data
structure (that is each face works on the cells surrounding it) the faces are
stored on one side of periodic boundary, but still have information about cells
on each side of periodicity.  It is enough to browse through one copy of the
periodic face for almost all numerical algorithms in T-Flows.

> **_Note:_** The exception is Lagrangian particle tracking which needs the
periodic face-pairs not to allow particles escape the computational domain.  
To visualize those period pairs, one can read the file ```lid_driven.shadows.vtu```.  
If plotted together with ```lid_driven.faces.vtu```, they close the domain.

### Running the simulation.  

In order to run a simulation, you have to compile _Process_ as explained in
[this section](#compiling_sub_programs).  In addition to that, you need a special
file called ```control``` which controls all the numerical parameters
(by _numerical_ we mean discretization), physical models (turbulence,
multiphase, multiple domains) and linear solver parameters (either for native
or PETSc, if the code was compiled with PETSc) and most importantly sets
boundary condition.  Although the amount of options which can be prescribed in
```control``` file is rather big, T-Flows' _Processs_ is very flexible and for
all the options which are not specified, it takes default values.  At its bare
minimum, ```control``` file needs problem name specified, to know which grids
to read (files ```.cfn``` and ```.dim```) and boundary conditions.

In order to facilitate the setup of ```control``` files, _Convert_ and _Generate_
create a control file template for the given grid.  In this case, the control
file template is called ```control_template_for_lid_driven``` and it resides in
the current directory.  Feel free to open it, because it comes with a lots of
useful information:
```
#-----------------------------------------------------------
# This is a minimum control file for T-Flows for the domain
# "lid_driven"
#
# One can use it as a staring point for defining entries
# in T-Flows' control file.  These files are created at
# the end of each call to Convert or Generate sub-programs.
#
# Each line starting with a # is, obviously, a comment.
# Lines starting with ! or % are also skipped as comments.
#-----------------------------------------------------------

#-----------------------------------------------------------
# Problem name must be specified and corresponds to the
# base name (without extensions) of the grid file used.
#-----------------------------------------------------------
 PROBLEM_NAME        lid_driven

#-----------------------------------------------------------
# Boundary conditions also must be specified and are
# listed here, for this grid, with some dummy values
#
# T-Flows neglects what it doesn't need to read, so in
# the lines below, only "wall" will be read as the type of
# boundary conditions, everything in () braces is skipped.
#
# The same holds for the values listed.  If simulation is
# laminar or LES, T-Flows will skip the values for "k",
# "eps", "zeta" and "f22"
#-----------------------------------------------------------
 BOUNDARY_CONDITION upper_wall
   TYPE             wall  (or: inflow / outflow / pressure / convective)
   VARIABLES        u     v     w     t    kin    eps    zeta     f22
   VALUES           0.0   0.0   0.0   10   1e-2   1e-3   6.6e-2   1e-3

 BOUNDARY_CONDITION side_walls
   TYPE             wall  (or: inflow / outflow / pressure / convective)
   VARIABLES        u     v     w     t    kin    eps    zeta     f22
   VALUES           0.0   0.0   0.0   10   1e-2   1e-3   6.6e-2   1e-3

 BOUNDARY_CONDITION lower_wall
   TYPE             wall  (or: inflow / outflow / pressure / convective)
   VARIABLES        u     v     w     t    kin    eps    zeta     f22
   VALUES           0.0   0.0   0.0   10   1e-2   1e-3   6.6e-2   1e-3

#-----------------------------------------------------------
# And, that's it, what is listed above is the bare minimum.
# If nothing else is prescribed, T-Flows will take default
# values as documented in Documents/all_control_keywords.
#
# Clearly, such a simulation wouldn't be very useful, so
# in the following are some of the most essential options
# for the control file. Simply uncomment what you need.
#
# All the entries (except the sub-entries such as: TYPE,
# VARIABLES/VALUES after each of the BOUNDARY_CONDITION
# can be inserted in any order.  T-Flows will take whatever
# it needs, whenever it needs it.  The entire boundary
# section may as well be at he end of the file.
#-----------------------------------------------------------

#-----------------------------------------------------------
# Time stepping
#-----------------------------------------------------------
! TIME_STEP                 0.001
! NUMBER_OF_TIME_STEPS   1200
! RESULTS_SAVE_INTERVAL   300  (how often will it save)

#-----------------------------------------------------------
# Physical properties
#-----------------------------------------------------------
! MASS_DENSITY           1.0
! THERMAL_CONDUCTIVITY   1.4e-4
! DYNAMIC_VISCOSITY      1.0e-4
! HEAT_CAPACITY          1.0

#-----------------------------------------------------------
# Monitoring points
#-----------------------------------------------------------
! NUMBER_OF_MONITORING_POINTS  4
!   MONITORING_POINT_001       0.0  0.0  2.0
!   MONITORING_POINT_002       0.0  0.0  2.5
!   MONITORING_POINT_003       0.0  0.0  3.0
!   MONITORING_POINT_004       0.0  0.0  3.5
! POINT_FOR_MONITORING_PLANES  0.17  0.17  1.51

#-----------------------------------------------------------
# Initial condition
#
# These are very similar to BOUNDARY_CONDITION section above
# and also here sub-entries VARIABLES and VALUES must be
# specified in the given order.
#-----------------------------------------------------------
! INITIAL_CONDITION
!   VARIABLES        u     v     w     t    kin    eps    zeta     f22
!   VALUES           0.0   0.0   0.0   10   1e-2   1e-3   6.6e-2   1e-3

#-----------------------------------------------------------
# For flows with inlets and outlets, initialization by a
# pressure-like potential field could prove to be useful
#-----------------------------------------------------------
!  POTENTIAL_INITIALIZATION      yes

#-----------------------------------------------------------
# For more customization, check the file:
# Documents/all_control_keywords.
#-----------------------------------------------------------
```
Please read through it, it gives a lots of explanations which are probably not
necessary to be repeated here agian.  We will cover different options in
subsequent sections of this manual.  For the sake of shortness, copy the
template control file to just ```control``` file with:
```
cp  control_template_for_lid_driven  control
```

remove all non-essential comments and set the velocity of the moving wall to
```1.0``` in the section ```BOUNDARY_CONDITION moving_wall``` until you get this:
```
# Problem name
 PROBLEM_NAME        lid_driven

# Boundary conditions
 BOUNDARY_CONDITION upper_wall
   TYPE             wall  
   VARIABLES        u     v     w
   VALUES           1.0   0.0   0.0

 BOUNDARY_CONDITION side_walls
   TYPE             wall
   VARIABLES        u     v     w
   VALUES           0.0   0.0   0.0

 BOUNDARY_CONDITION lower_wall
   TYPE             wall  
   VARIABLES        u     v     w
   VALUES           0.0   0.0   0.0
```

At this point, you are ready to run.  Invoke _Process_ by issuing command:
```
../../../../Binaries/Process > out  &
```

Since _Process_ writes a lot of information on the screen while it is computing,
it is useful to re-direct the output to a log file, here simply called ```out```.  
We also send the process in the background with an ampersand ```&``` at the end
of the command line.  Next, let's analyze the output from _Process_.  
It starts with a header:
```
 #=======================================================================
 #
 #    ______________________ ____    ________  __      __  _________
 #    \__    ___/\_   _____/|    |   \_____  \/  \    /  \/   _____/
 #      |    |    |    __)  |    |    /   |   \   \/\/   /\_____  \
 #      |    |    |     \   |    |___/    |    \        / /        \
 #      |____|    \___  /   |_______ \_______  /\__/\  / /_______  /
 #                    \/            \/       \/      \/          \/
 #                      ___
 #                     / _ \_______  _______ ___ ___
 #                    / ___/ __/ _ \/ __/ -_|_-<(_-<
 #                   /_/  /_/  \___/\__/\__/___/___/
 #
 #                        Double precision mode
 #                     Compiled with PETSc solvers
 #-----------------------------------------------------------------------
```

which gives two important pieces of information at the end: that the _Process_
was compiled in 64-bit representation of floating point numbers (```Double
precision mode```) and that it was linked with PETSc libraries (```Compiled
with PETSc solvers```).  What follows is a bunch of parameters it reads from
the ```control``` file, and it is generally worth observing the lines such as
these (some are omitted):
```
 # NOTE! Could not find the keyword: NUMBER_OF_TIME_STEPS. Using the default:      1200
 # ...
 # NOTE! Could not find the keyword: HEAT_TRANSFER. Using the default: no        
 # ...
 # NOTE! Could not find the keyword: TURBULENCE_MODEL. Using the default: none   
 # NOTE! Could not find the keyword: INTERFACE_TRACKING. Using the default: no    
 # NOTE! Could not find the keyword: PARTICLE_TRACKING. Using the default: no    
 # NOTE! Could not find the keyword: TIME_STEP. Using the default: 1.000E-02     
```

At this point, most important is that, because neither the number of time steps i
nor the time step is prescribed in ```control``` file, it sets them to default
values of ```1200``` and ```1.0E-2``` respectively.

Each time step performed by _Process_ prints some information on convergence of
solution procedure, and it looks as follows:
```
                                         #===============================================#
                                         #                                               #
                                         #              Time step :    131               #
                                         #                                               #
                                         #    Simulation time : 1.310E+00 [s]            #
                                         #    Wall-clock time : 000:00:02 [hhh:mm:ss]    #
                                         #                                               #
                                         #-----------------------------------------------#

  #=============================================================================================================================#
  #                                                        Iteration:  1                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  2 2.966E-08 | V   :  0 0.000E+00 | W   :  2 4.583E-09 | PP  : 14 7.786E-07 | MASS:    1.476E-08 |                    #
  #=============================================================================================================================#
  #                                                        Iteration:  2                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  2 1.390E-08 | V   :  0 0.000E+00 | W   :  2 2.981E-09 | PP  : 14 8.401E-07 | MASS:    1.706E-08 |                    #
  #=============================================================================================================================#
  #                                                        Iteration:  3                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  2 6.594E-09 | V   :  0 0.000E+00 | W   :  1 7.434E-07 | PP  : 14 7.609E-07 | MASS:    1.528E-08 |                    #
  #--------------------+=========================================+=========================================+--------------------#
                       #    Maximum Courant number: 1.602E-01    |    Maximum Peclet number: 4.006E+00     #
                       #---------------------------+-------------+-------------+---------------------------#
                       #    Flux x :  0.000E+00    |    Flux y :  -1.49E-10    |    Flux z :   0.00E+00    #
                       #    Pdrop x:  0.000E+00    |    Pdrop y:   0.00E+00    |    Pdrop z:   0.00E+00    #
                       #---------------------------+---------------------------+---------------------------#
```
In the header, it shows the current time step, current simulation (numerical)
time and elapsed wall-clock time.  That is followed by information on achieved
convergence within each time step.  The equations being solved for this case
are three velocity components (```U```, ```V``` and ```W```) and pressure
correction (```PP```).  For each of them it shows number of iterations performed
in linear solver (the small integer values), as well as the final residual
achieved (the number written in exponential form).  The footer shows the
maximum Courant and Peclet (_Pe_) number at the end of the time step, as well
as global fluxes and pressure drops in the current time step.  For this
simulation, global fluxes and pressure drops are of no relevance, but maximum
Courant number gives us insight about convergence.  If it grows beyond a hundred,
we will likely face divergence.  

Since we know that this test case reaches a steady solution in the end, we should
expect number of iterations to drop to zero.  In this case, it is not quite like
that.  Last performed time step shows:
```
                                         #===============================================#
                                         #                                               #
                                         #              Time step :   1200               #
                                         #                                               #
                                         #    Simulation time : 1.200E+01 [s]            #
                                         #    Wall-clock time : 000:00:19 [hhh:mm:ss]    #
                                         #                                               #
                                         #-----------------------------------------------#

  #=============================================================================================================================#
  #                                                        Iteration:  1                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  1 2.462E-08 | V   :  0 0.000E+00 | W   :  1 1.886E-08 | PP  :  2 5.501E-07 | MASS:    1.053E-08 |                    #
  #=============================================================================================================================#
  #                                                        Iteration:  2                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  1 1.134E-08 | V   :  0 0.000E+00 | W   :  1 8.567E-09 | PP  :  0 0.000E+00 | MASS:    1.716E-08 |                    #
  #=============================================================================================================================#
  #                                                        Iteration:  3                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  1 5.195E-09 | V   :  0 0.000E+00 | W   :  1 3.848E-09 | PP  :  2 5.312E-07 | MASS:    1.184E-08 |                    #
  #--------------------+=========================================+=========================================+--------------------#
                       #    Maximum Courant number: 1.672E-01    |    Maximum Peclet number: 4.179E+00     #
                       #---------------------------+-------------+-------------+---------------------------#
                       #    Flux x :  0.000E+00    |    Flux y :  -2.07E-11    |    Flux z :   0.00E+00    #
                       #    Pdrop x:  0.000E+00    |    Pdrop y:   0.00E+00    |    Pdrop z:   0.00E+00    #
                       #---------------------------+---------------------------+---------------------------#
```
Number of iterations are all one or two, but not quite zero yet.  So, if we
want to be purists, we can't say we achieved steady state to the tolerance set
by T-Flows.  In order to work around it, we could either increase the number
of time steps, or the time step itself.  Since CFL is well below unity, let's
re-run the simulation with ```control``` file modified as:
```

# Problem name
 PROBLEM_NAME        lid_driven

# Time step
 TIME_STEP  0.1

# Boundary conditions
 BOUNDARY_CONDITION upper_wall
   TYPE             wall
   VARIABLES        u     v     w
   VALUES           1.0   0.0   0.0

 BOUNDARY_CONDITION side_walls
   TYPE             wall
   VARIABLES        u     v     w
   VALUES           0.0   0.0   0.0

 BOUNDARY_CONDITION lower_wall
   TYPE             wall
   VARIABLES        u     v     w
   VALUES           0.0   0.0   0.0
```

With time step set to ```0.1```, convergence is observed after 312 time steps
since all iterations for all computed quantities fell to zero:
```
                                         #===============================================#
                                         #                                               #
                                         #              Time step :    312               #
                                         #                                               #
                                         #    Simulation time : 3.120E+01 [s]            #
                                         #    Wall-clock time : 000:00:05 [hhh:mm:ss]    #
                                         #                                               #
                                         #-----------------------------------------------#

  #=============================================================================================================================#
  #                                                        Iteration:  1                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  0 9.789E-07 | V   :  0 0.000E+00 | W   :  0 9.640E-07 | PP  :  0 0.000E+00 | MASS:    1.856E-07 |                    #
  #=============================================================================================================================#
  #                                                        Iteration:  2                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  0 9.789E-07 | V   :  0 0.000E+00 | W   :  0 9.640E-07 | PP  :  0 0.000E+00 | MASS:    1.856E-07 |                    #
  #=============================================================================================================================#
  #                                                        Iteration:  3                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  0 9.789E-07 | V   :  0 0.000E+00 | W   :  0 9.640E-07 | PP  :  0 0.000E+00 | MASS:    1.856E-07 |                    #
  #--------------------+=========================================+=========================================+--------------------#
                       #    Maximum Courant number: 1.672E+00    |    Maximum Peclet number: 4.181E+00     #
                       #---------------------------+-------------+-------------+---------------------------#
                       #    Flux x :  0.000E+00    |    Flux y :   8.54E-11    |    Flux z :   0.00E+00    #
                       #    Pdrop x:  0.000E+00    |    Pdrop y:   0.00E+00    |    Pdrop z:   0.00E+00    #
                       #---------------------------+---------------------------+---------------------------#
```
That is it, you got your first CFD solutions with T-Flows.  You can read the
final results (```lid_driven-ts001200.vtu```) in ParaView, and explore different
options for visualization of results.  We chose to plot a cut-plane through the
middle of the domain vectors represented glyphs and pressure with a colormap:

![Lid-driven hexa solution!](Documentation/Manual/Figures/lid_driven_hexa_solution.png "Lid-driven hexa solution")

but you obviously have the freedom to explore other options in ParaView.

#### Things to try next

##### Refining the grid

This grid was rather coarse, only 20 x 20 cells in the cross plane, and three
rows of cells in the spanwise direction.  Since spanwise direction is periodic,
it doesn't make sense in it, but you can refine in other directions.  To that
end, you should modify the ```lid_driven.geo``` file in line which currently reads:
```
NA = 20;    // resolution in length and height
```
to:
```
NA = 40;    // resolution in length and height
```
After that, you should:
- Re-run GMSH with: ```gmsh -3 lid_driven.geo```
- Re-run _Convert_ with the newly created ```lid_driven.msh``` file
- Re-Run _Process_ with ```../../../../Binaries/Process > out_finer_mesh  &```

and see how the results look on a refined mesh.  You can carry on my increasing
the parameter ```NA``` to 80, 160, or even higher values.

#### Typical workflow

The sequence outlined above represents a typical workflow in T-Flows, hence:
1. Set or modify parameters in ```.geo``` file
2. Re-run GMSH
3. Re-run _Convert_ to get the files in the T-Flows format
4. Re-run the simulation
5. Check results, go back to step 1, if further improvements are needed.

#### Advice on speeding up the conversion

Since in the workflow outline above, entering values into _Convert_ every
time you run can become tedious, you can speed up the process by creating a
special file, we call it ```convert.scr``` which contains all the answers you
give to _Convert_.  In this particular case, the file ```convert.scr```
would look like:
```
lid_driven.msh
no
skip
2
4
skip
skip
skip
```

Once you have this file, you can invoke _Convert_ with:
```
../../../../Binaries/Convert < convert.scr
```
and don't have to type anything.

#### A utility to shorten those paths

The T-Flows directory structure, test cases are many and are categorized in
many directories and their directories.  Invoking executables such as
_Convert_ and _Process_ can be tedious and error prone with all those dots and
slashes before the executable name.  In order to avoid it, directory ```[root]/Sources/Utilites```
holds a bash script called ```seek_binaries.sh```.  It is an executable script
and if you place it somewhere in your path (say in ```/usr/local/bin```) every
time you call it from a test directory, it will seek the executables and make
soft links in your test directory.  For this particular case it will create links:
```
Hexa/
├── Convert -> ../../../../Binaries/Convert
└── Process -> ../../../../Binaries/Process
```
and from there on you can invoke them with simple ```./Convert``` and
```./Process```.

#### Try some of the standard T-Flows's cases

Better and more elaborate test cases for the lid-driven cavity flow have been
set in ```[root]/Tests/Laminar/Cavity/Lid_Driven```.  Feel free to explore them further.

### On polyhedral grid <a name="demo_lid_driven_dual"></a>

During the conversion process outlined above, you were asked if you wanted to
created a dual grid which we simply skipped.  In this section, we will show you
what is behind it.  In order to run this case, please go to the directory
```[root]/Tests/Manual/Lid_Driven_Cavity/Dual``` where you will find the
following files:
```
Dual/
├── convert.scr
├── lid_driven.geo
└── lid_driven.msh.gz
```

> **_Note:_** The file ```lid_driven.msh.gz``` is here for the same reason as
explained above: to keep you going even if you don't have GMSH installed.

The ```.geo``` file is almost the same as the one we used for orthogonal grid.  
If you run the two ```.geo``` files through ```diff``` command:
```
diff lid_driven.geo ../Hexa/lid_driven.geo
```
you will see the following output:
```
> // Set surface to be transfinite and quadrilateral
> Transfinite Surface {1} = {1, 2, 3, 4};  Recombine Surface {1};
```
meaning that only the line which instructs GMSH to create quadrilateral mesh
is missing.  Since that is the only change and boundary conditions are exactly
the same, the ```convert.scr``` is the same as in previous case, we can re-use it.

Without being instructed to create quadrilateral cells, GMSH creates triangular.  
If you run this script:
```
gmsh -3 lid_driven.geo
```
and convert it to T-Flows format with:
```
./Convert < convert.scr
```
it will creates a computation grid like this:

![Lid-driven prismatic!](Documentation/Manual/Figures/lid_driven.png "Lid-driven prismatic grid")

For a node-base numerical framework this grid looks rather descent, but T-Flows
is cell-based and on grids based on triangular prisms (and tertahedra) induce:
- poor accuracy of gradient computation which impacts the accuracy of the overall
  numerical scheme, and
- large explicit diffusion terms causing slower convergence of pressure-velocity
  coupling algorithm (SIMPLE or PISO in T-Flows)

The disadvantage of the tetrahedral grids has been recognized long ago, and
polyhedral grids have been introduced as their alternative.  In mathematical
sense, a polyhedral is nothing more than a [dual graph](https://en.wikipedia.org/wiki/Dual_graph) of
tetrahedral grid, and that's why _Convert_ calls this a _dual_ grid.  

Having mentioned the duality of the two grids, it is worth briefly explaining
how _Convert_ performs this.  It reads a grid with triangular prisms and/or
tetrahedra in the usual way, does its conversion (finds the connectivity between
cells, faces, nodes and edges it needs) and in the next step creates a graph
dual of the first mesh.  The only differnce for you, as a user is that when
prompted by the question:
```
 #=================================================
 # Would you like to create a dual grid? (yes/no)
 #-------------------------------------------------
```
you answer ```yes```.  From that point on, everything is the same as explained
in the [section above](#demo_lid_driven_hexa_convert).  To save you from
typing all the answers, the file ```convert.scr``` is also provided in the
current directory.  If you created such a file for the case with
[hexaderal cells](#demo_lid_driven_hexa), it would be almost the same as
this one.  If you ran two files through a diff command:
```
diff convert.scr ../Hexa/convert.scr
```
you would see only this:
```
< yes
---
> no
```
the only difference between the files is the fact that you answered ```yes```
in this case when aasked whether you want to create a dual grid.

Anyhow, in order to distinguish between the original and the dual grid,
_Convert_ adds extension ```_dual``` to all the file names it creates.  
Thus, after running the _Convert_ you will have the following files in the directory:
```
Dual/
├── lid_driven_dual.cfn
├── lid_driven_dual.dim
├── lid_driven_dual.edges.vtu
├── lid_driven_dual.faces.vtu
├── lid_driven_dual.shadows.vtu
└── lid_driven_dual.vtu
```

The first two are T-Flows's native files for further processing, and the
remaining four are for you to explore with Paraview.  Opening the
```lid_driven_dual.shadows.vtu```, shows boundary conditions, which are the
same as in the [previous case](#demo_lid_driven_hexa):

![Lid-driven polyhedral!](Documentation/Manual/Figures/lid_driven_dual.png "Lid-driven polyhedral grid")

To run the case, we have already provided the ```control``` file, derived from
the previous case we ran.  Here it is in full:
```
# Problem
PROBLEM_NAME        lid_driven_dual

# Time stepping
TIME_STEP  0.1

# Boundary conditions
BOUNDARY_CONDITION upper_wall
  TYPE             wall
  VARIABLES        u     v     w
  VALUES           1.0   0.0   0.0

BOUNDARY_CONDITION side_walls
  TYPE             wall
  VARIABLES        u     v     w
  VALUES           0.0   0.0   0.0

BOUNDARY_CONDITION lower_wall
  TYPE             wall
  VARIABLES        u     v     w
  VALUES           0.0   0.0   0.0
```
The only novelty compared to the previous case is line with the ```PROBLEM_NAME```.  
It is set to ```lid_driven_dual``` here.  In any case, after running the
simulation, a possible representation fo the solution looks like:

![Lid-driven hexa solution!](Documentation/Manual/Figures/lid_driven_dual_solution.png "Lid-driven dual solution")

## Thermally-driven cavity flow <a name="demo_thermally_driven"></a>

Thermally-driven cavity flow bears many similarities with the lid-driven
[cavity flow](#tests_lid_driven).  Both of these flows occur in enclosures with
square cross-section, both are without inflows and outflows facilitating
prescription of boundary conditions, both are occuring in ddomains which are
long enough in spanwise direction so that the assumption of two-dimsionality or
periodicity can be made.  Owing to their simplicity, both of these cases have
widely been used by CFD community for benchmarking and verification of CFD
codes and both are well documented.

The biggest differences between the cases is the driving force.  Whereas the
lid-driven cavity flow is driven by shear created by top moving wall,
thermally-driven cavity is driven by the buoyancy forces occurring on vertical
opposing sides of the problem domain:

![Thermally-driven cavity domain!](Documentation/Manual/Figures/thermally_driven_domain.png "Thermally-driven cavity domain")

In the figure above the left (red) wall is kept at higher temperature than the
right wall (blue), which creates clockwise motion of the fluid.  

We will cover the thermall-driven cavity flow in two modelling ways: with
Boussinesq approximation, and with variable physical properties.

### With Boussinesq approximation <a name="demo_thermally_driven_boussinesq"> </a>

We will use this case to introduce a few new concepts in T-Flows:
- instruct _Process_ to solve energy equation (temperature)
- specify different boundary conditions for temperature (Dirichlet or Neumann)
- changing from default under-relaxation parameters
- setting the linear solver tolerances
- setting the initial condition
To solve this case, please go to the directory: ```[root]/Tests/Manual/Thermally_Driven/Direct```
where you can find the following files:

```
Direct/
├── convert.scr
├── therm_driven.geo
└── therm_driven.msh.gz
```

The ```.geo``` file is the GMSH script.  Since the geometries for the lid-driven
[cavity flow](#tests_lid_driven) and this case are almost the same, the
similarity in the ```.geo``` shouldn't be a suprise.  You can find the
following differences:
```
A    =  1.0;  // length and height of the cavity
B    =  0.1;  // width of the cavity
NA   = 60;    // resolution in length and height
NB   =  3;    // resolution in width (periodic direction)
BUMP = 0.1;   // control clustering towards the walls
...
...
// Set lines to transfinite
// (+1 is because GMSH expects number of points here)
Transfinite Curve {1, 2, 3, 4} = NA+1  Using Bump BUMP;
...
...
// Define boundary conditions
Physical Surface("upper_wall", 27) = {21};
Physical Surface("left_wall", 28) = {17};
Physical Surface("right_wall", 29) = {25};
Physical Surface("lower_wall", 30) = {13};
Physical Surface("periodic_y", 31) = {1, 26};
...
...
```

The novelties include:
- increased resolution ```NA```
- new parameter ```BUMP``` which ...
- controls the clustering of the lines towards the walls in the
- call to ```Transfinite Curve```
- boundary condition names have changed; ```left_wall``` and ```right_wall```
  are new, the ```side_walls``` has been dropped.

First thing would be to generate the mesh with:
```
gmsh -3 therm_driven.geo
```

which will create ```therm_driven.msh```.

> **_Note:_** For those without GMSH, the file ```therm_driven.msh.gz``` is provided.

We advise to run ```seek_binaries.sh``` next to get all executables to current
folder, after which you can call _Convert_ with provided script:
```
./Convert < convert.scr
```

> **_Note:_** The script ```convert.scr``` is also almost the same as in the
case of the lid-driven cavity.  You can explore it yourself to check that only
the input file name is different (```therm_driven.msh``` instead of
```lid_driven.msh```), and that periodic boundary is listed under different
number (```5```, not ```4```).

If you visualize the file ```therm_driven.faces.vtu``` created by _Convert_ you
will see this:

![Thermally-driven front!](Documentation/Manual/Figures/thermally_driven_front.png "Thermally-driven front")

a grid with cells clustered towards the walls for better resolution of boundary
layers.  It is not a waste of time to visualize grids created by _Convert_ to
make sure that no anomalities are present.  

With grids converted, we should set up the ```control``` file.  To that end,
lets' consider the section related to boundary conditions created by _Convert_
(```control_template_for_therm_driven```, with all variables which are not solved
removed for the sake of clarity:
```
BOUNDARY_CONDITION upper_wall
  TYPE             wall
  VARIABLES        u     v     w     q
  VALUES           0.0   0.0   0.0   0.0

BOUNDARY_CONDITION left_walls
  TYPE             wall
  VARIABLES        u     v     w     t
  VALUES           0.0   0.0   0.0   1.0

BOUNDARY_CONDITION right_walls
  TYPE             wall
  VARIABLES        u     v     w     t
  VALUES           0.0   0.0   0.0   0.0

BOUNDARY_CONDITION lower_wall
  TYPE             wall
  VARIABLES        u     v     w     q
  VALUES           0.0   0.0   0.0   0.0

```
In the above, we left velocities at all walls to zero, we set temperature at
the left wall to ```1.0```, temperature at the ```right_wall``` to ```0.0```,
but in order to specify that ```upper_wall``` and ```lower_wall``` are with
prescribed heat flux rather than temperature, we change letter ```t``` to ```q```
(as a usual symbol for heat flux).  We also set it to zero, because they are insulated.

We should also instruct _Process_ to solve for temperature, which is obtained with the line:
```
HEAT_TRANSFER    yes
```
Given that the temperatures at the boundaries are zero and one, and since there
are no internal sources or sinks, we may expect final average temperature to be
around 0.5.  Hence, it makes sense to set initial temperature to be 0.5 everywhere:
```
 INITIAL_CONDITION
   VARIABLES           u     v     w     t
   VALUES              0.0   0.0   0.0   0.5
```

Next, we should also specify physical properties.  For this case, we will solve
equations in their non-dimensional form:

EQUATIONS WHICH ARE MISSING

from which we can see that the flow is fully characterized with two
non-dimenensional numbers: Rayleigh (_Ra_) and Prandtl (_Pr_).  _Process_
doesn't know about these numbers, so we set physical properties in a way to
reach the desired values.  If we want to solve the case for:
- _Ra_ = 1.0e+6
- _Pr_ = 0.71

physical properties section in the control file should read:
```
 # Properties based on Pr and Ra numbers:
 # Pr = 0.71
 # Ra = 10e6
 # mu     = 1.0 / sqrt(Pr * Ra) = 0.001186781658
 # lambda = sqrt(Pr / Ra)       = 0.000842614977
 DYNAMIC_VISCOSITY      0.001186781658
 THERMAL_CONDUCTIVITY   0.000842614977
```

Mass density, heat capacity and thermal expansion coefficient are not set, and
_Process_ will set them to their default values of 1.0.

> **_Note:_** Be reminded that default values for all parameters needed by
_Process_ are outlined in the file: ```[root]/Documentation/all_control_keywords```.

We should also instruct _Process_ that we want to use Boussinesq approximation
to solve the system, which is obtained with the line:
```
 BUOYANCY         thermal
```

_Process_ also needs to know the direction and magnitude of the gravitational
vector which is set by:
```
GRAVITATIONAL_VECTOR   0.0  0.0  -1.0
```

Since we know in advance that the velocities in this case are rather small we
can increase the time step to 1.0:
```
TIME_STEP    1.0
```
and, since we seek a steady solution, we can reduce the saving frequency as:
```
RESULTS_SAVE_INTERVAL    300
```

_Process_, by default, links velocities and pressure through the SIMPLE algorith.  
(The other option in the code is PISO.)  An important aspect of the
pressure-velocity algorithms are the _under-relaxation_ factors, which may be
fiddled with in order to improve convergence of the solution procedure.  If not
specified, _Process_ sets very conservative (read: small) values to make sure
the system will converge.  For this case, given that it's relativelly simple,
laminar, and that we have to make sure that temperatures reach final
stratification, we may want to set them explicitly.  From simulations we
conducted beforehand, we found these values to work for this case:
```
 PRESSURE_MOMENTUM_COUPLING             simple
 SIMPLE_UNDERRELAXATION_FOR_MOMENTUM    0.7
 SIMPLE_UNDERRELAXATION_FOR_ENERGY      0.3
 SIMPLE_UNDERRELAXATION_FOR_PRESSURE    0.8
```

> **_Note:_** The first line is optional, since default pressure-velocity
coupling in T-Flows is SIMPLE.

Another thing worth noting for this case is that the default linear solver
parameters might not be the best ones (for all variables solved, it is 1.0e-6).  
This may be too tight for velocities and temperature, and a bit too loose for
pressure.  We therefore set them as following:
```
 LINEAR_SOLVERS                     native
 TOLERANCE_FOR_MOMENTUM_SOLVER      1.e-3
 TOLERANCE_FOR_ENERGY_SOLVER        1.e-3
 TOLERANCE_FOR_PRESSURE_SOLVER      1.e-7
```
First line tells _Process_ to use its owen (_native_) solvers, whereas the
remaining three lines are self-explanatory, we trully believe.

With all this in place, the entire control file may read like this:
```
# Problem name
 PROBLEM_NAME             therm_driven

# Related to temperature
 HEAT_TRANSFER            yes
 BUOYANCY                 thermal
 GRAVITATIONAL_VECTOR     0.0  0.0  -1.0

# Time stepping
 TIME_STEP                1
 RESULTS_SAVE_INTERVAL  300

# Properties based on Pr and Ra numbers:
# Pr = 0.71
# Ra = 10e6
# mu     = 1.0 / sqrt(Pr * Ra) = 0.001186781658
# lambda = sqrt(Pr / Ra)       = 0.000842614977
 DYNAMIC_VISCOSITY      0.001186781658
 THERMAL_CONDUCTIVITY   0.000842614977

# Initial conditions
 INITIAL_CONDITION
   VARIABLES           u     v     w     t
   VALUES              0.0   0.0   0.0   0.5

# Numerical parameters
 PRESSURE_MOMENTUM_COUPLING             simple
 SIMPLE_UNDERRELAXATION_FOR_MOMENTUM    0.7
 SIMPLE_UNDERRELAXATION_FOR_ENERGY      0.3
 SIMPLE_UNDERRELAXATION_FOR_PRESSURE    0.8
 TOLERANCE_FOR_SIMPLE_ALGORITHM         1.0e-3

# Linear solver parameters
 LINEAR_SOLVERS                     native
 TOLERANCE_FOR_MOMENTUM_SOLVER      1.e-3
 TOLERANCE_FOR_ENERGY_SOLVER        1.e-3
 TOLERANCE_FOR_PRESSURE_SOLVER      1.e-7

# Initial condition
 INITIAL_CONDITION
   VARIABLES           u     v     w     t
   VALUES              0.0   0.0   0.0   0.5

# Boundary conditions
 BOUNDARY_CONDITION upper_wall
   TYPE             wall
   VARIABLES        u     v     w     q
   VALUES           0.0   0.0   0.0   0.0

 BOUNDARY_CONDITION left_wall
   TYPE             wall
   VARIABLES        u     v     w     t
   VALUES           0.0   0.0   0.0   1.0

 BOUNDARY_CONDITION right_wall
   TYPE             wall
   VARIABLES        u     v     w     t
   VALUES           0.0   0.0   0.0   0.0

 BOUNDARY_CONDITION lower_wall
   TYPE             wall
   VARIABLES        u     v     w     q
   VALUES           0.0   0.0   0.0   0.0
```
We say "may" because, as mentioned above, the order of individual entries
doesn't matter.  Feel free to copy the above contents to ```control``` file and
launch the simulation with:
```
./Process > out &
```
given that you created soft links to executables in this directory with
```seek_binaries.sh``` script.

With these settings in the control file, you will reach convergence in about
1000 time step, as obvious from the output:
```
                                         #===============================================#
                                         #                                               #
                                         #              Time step :   1000               #
                                         #                                               #
                                         #    Simulation time : 1.000E+03 [s]            #
                                         #    Wall-clock time : 000:01:39 [hhh:mm:ss]    #
                                         #                                               #
                                         #-----------------------------------------------#

  #=============================================================================================================================#
  #                                                        Iteration:  1                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  0 9.044E-04 | V   :  0 0.000E+00 | W   :  0 8.232E-04 | PP  :  0 0.000E+00 | MASS:    9.222E-07 | T   :  0 9.996E-04 #
  #=============================================================================================================================#
  #                                                        Iteration:  2                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  0 9.044E-04 | V   :  0 0.000E+00 | W   :  0 8.232E-04 | PP  :  0 0.000E+00 | MASS:    9.222E-07 | T   :  0 9.996E-04 #
  #=============================================================================================================================#
  #                                                        Iteration:  3                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  0 9.044E-04 | V   :  0 0.000E+00 | W   :  0 8.232E-04 | PP  :  0 0.000E+00 | MASS:    9.222E-07 | T   :  0 9.996E-04 #
  #--------------------+=========================================+=========================================+--------------------#
                       #    Maximum Courant number: 7.483E+00    |    Maximum Peclet number: 5.288E+00     #
                       #---------------------------+-------------+-------------+---------------------------#
                       #    Flux x :  0.000E+00    |    Flux y :   1.91E-12    |    Flux z :   0.00E+00    #
                       #    Pdrop x:  0.000E+00    |    Pdrop y:   0.00E+00    |    Pdrop z:   0.00E+00    #
                       #---------------------------+---------------------------+---------------------------#
```
Having obtained steady solution for this case, you may visualise some results
in ParaView by opening the file: ```therm_driven-ts001200.vtu```.  Here we show
solutions for temperature with velocity vectors scaled by their magitude:

![Thermally-driven solution!](Documentation/Manual/Figures/thermally_driven_solution.png "Thermally-driven solution")

#### Thing to try next

Better and more elaborate test cases for the thermall-driven cavity flow have
been placed in ```[root]/Tests/Laminar/Cavity/Thermally_Driven/Direct```, for a
range of _Ra_ numbers.  Feel free to explore them further.

### With variable physical properties <a name="demo_thermally_driven_variable"> </a>

Using Boussinesq hypotehsis is not the only way we can deal with buoyancy driven flows.
The alternative would be to be to change air density as the function of temperature, and
imposing a gravitational vector, _Process_ would work out buoyancy forces acting on
momentum equations.  Dependency of density on temperature has to be imposed in some way
and we see it as a look-up table.  We would have to use a _user function_ at each time
step to update the values of physical properties.

With this case we want to demonstrate:
- how to run a buoyancy driven flows without Boussinesq approximation
- how to write user functions.

A case which demonstrates how it is done resides in
```[root]/Tests/Manual/Thermally_Driven/Variable_Properties```.
The directory contains the following files:
```
[root]/Tests/Manual/Thermally_Driven/Variable_Properties
├── air.geo
├── air_properties_at_1_bar.dat
├── control
├── convert.scr
├── readme
└── User_Mod/
    ├── Beginning_Of_Simulation.f90
    ├── Beginning_Of_Time_Step.f90
    └── Types.f90
```
The purpose of files ```air.geo```, ```control``` and ```convert.scr```
should be clear by now.  If not, refer to previous sections on lid- and
thermally-driven cavities.
What is new is the look-up table in file ```air_properties_at_1_bar.dat```.  
If you open it, you can see that it lists air properties for a range of
temperatures from -20 to 125 in increments of 5 degrees.  

In addition to these, there is a sub-directory called ```User_Mod``` with three
functions.  The way T-Flows handles user functions is as follows.  In the
folder which holds sources for _Process_, there is also a sub-folder called
```User_Mod``` holding a number of functions user can modify.  They are:
```
[root]/Sources/Process/User_Mod
├── Before_Exit.f90
├── Beginning_Of_Compute_Energy.f90
├── Beginning_Of_Compute_Momentum.f90
├── Beginning_Of_Compute_Pressure.f90
├── Beginning_Of_Compute_Scalar.f90
├── Beginning_Of_Compute_Vof.f90
├── Beginning_Of_Correct_Velocity.f90
├── Beginning_Of_Iteration.f90
├── Beginning_Of_Simulation.f90
├── Beginning_Of_Time_Step.f90
├── Calculate_Mean.f90
├── End_Of_Compute_Energy.f90
├── End_Of_Compute_Momentum.f90
├── End_Of_Compute_Pressure.f90
├── End_Of_Compute_Scalar.f90
├── End_Of_Compute_Vof.f90
├── End_Of_Correct_Velocity.f90
├── End_Of_Iteration.f90
├── End_Of_Simulation.f90
├── End_Of_Time_Step.f90
├── Force.f90
├── Initialize_Variables.f90
├── Insert_Particles.f90
├── Interface_Exchange.f90
├── Save_Results.f90
├── Save_Swarm.f90
├── Source.f90
└── Types.f90
```
Their names are self-explanatory, as their names correspond to places in
_Processs_ from which they are called.  When _Process_ is compiled with simple
```make``` command, the functions from the ```[root]/Sources/Process/User_Mod```
will be called and they, in essence, do nothing.  They are _empty hooks_ so to
speak.  If _Process_ is compiled with ```DIR_CASE=full_or_relative_path_to_case```,
sources (_empty hooks_) from ```[root]/Sources/Process/User_Mod``` will be
replaced by the links to ```User_Mod``` directory residing with the case.

Try to see it in action.  Change the current directory to
```[root]/Sources/Process``` and issue the commands:
```
make clean
make DIR_CASE=../../Tests/Manual/Thermally_Driven/Varible/
```
> **_Tip:_** When dealing with user functions, it is not a bad idea to conduct
```make clean``` every once in a while.  It is possible, and often the case,
that objects from the last compilation of the _Process_ are newer than user
function specified in the case.

When the compilation completes, the contents of the ```User_Mod``` in
```[root]/Sources/Process``` will read:
```
[root]/Sources/Process/User_Mod
├── Before_Exit.f90
├── Beginning_Of_Compute_Energy.f90
├── Beginning_Of_Compute_Momentum.f90
├── Beginning_Of_Compute_Pressure.f90
├── Beginning_Of_Compute_Scalar.f90
├── Beginning_Of_Compute_Vof.f90
├── Beginning_Of_Correct_Velocity.f90
├── Beginning_Of_Iteration.f90
├── Beginning_Of_Simulation.f90 -> ../../../Tests/Manual/Thermally_Driven/Varible/User_Mod/Beginning_Of_Simulation.f90
├── Beginning_Of_Time_Step.f90 -> ../../../Tests/Manual/Thermally_Driven/Varible/User_Mod/Beginning_Of_Time_Step.f90
├── Calculate_Mean.f90
├── End_Of_Compute_Energy.f90
├── End_Of_Compute_Momentum.f90
├── End_Of_Compute_Pressure.f90
├── End_Of_Compute_Scalar.f90
├── End_Of_Compute_Vof.f90
├── End_Of_Correct_Velocity.f90
├── End_Of_Iteration.f90
├── End_Of_Simulation.f90
├── End_Of_Time_Step.f90
├── Force.f90
├── Initialize_Variables.f90
├── Insert_Particles.f90
├── Interface_Exchange.f90
├── Save_Results.f90
├── Save_Swarm.f90
├── Source.f90
└── Types.f90 -> ../../../Tests/Manual/Thermally_Driven/Varible/User_Mod/Types.f90
```
The three files which are defined in the case's ```User_Mod``` are now linked
to the _Process'_ ```User_Mod```.  That means they are not _empty hooks_ any
longer, but whatever user modified them to do.

Having described the mechanism by which _Process_ deals with user functions,
we can describe each of them more closely.

The source ```Types.f90``` holds definion of new types which might be used with
user functions.  It is not always necessary to define new types for user
functions, but in this particular case we want to store look-up tables with variable air properties in memory.  It reads:
```
  1 !==============================================================================!
  2 !   Introduce new types to be used with User_Mod                               !
  3 !==============================================================================!
  4
  5   integer, parameter :: N_ITEMS = 30
  6
  7   real :: air_t     (N_ITEMS)
  8   real :: air_rho   (N_ITEMS)
  9   real :: air_mu    (N_ITEMS)
 10   real :: air_cp    (N_ITEMS)
 11   real :: air_lambda(N_ITEMS)

```
We know that there are 30 entries in the file ```air_properties_at_a_bar.dat```
and we therefore introduce a parameter called ```N_ITEMS``` in line 5 and set
it to ```30```.  Arrays holding phyiscal properties are statically allocate
with that size in lines 7 - 11.

The source ```Beginning_Of_Simulation``` is, clearly enough, called at the
beginning of simulation, and we wrote in a way to read physical properties:
```
  1 !==============================================================================!
  2   subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm, n, time)
  3 !------------------------------------------------------------------------------!
  4 !   This function is called at the beginning of simulation.                    !
  5 !------------------------------------------------------------------------------!
  6   implicit none
  7 !---------------------------------[Arguments]----------------------------------!
  8   type(Field_Type),    target :: Flow
  9   type(Turb_Type),     target :: Turb
 10   type(Vof_Type),      target :: Vof
 11   type(Swarm_Type),    target :: Swarm
 12   integer, intent(in)         :: n     ! time step
 13   real,    intent(in)         :: time  ! physical time
 14 !-----------------------------------[Locals]-----------------------------------!
 15   type(Grid_Type), pointer :: Grid
 16   integer                  :: i, fu
 17 !==============================================================================!
 18
 19   ! Take aliases
 20   Grid => Flow % pnt_grid
 21
 22   !----------------------------------------!
 23   !   Open file with physical properties   !
 24   !----------------------------------------!
 25   call File % Open_For_Reading_Ascii("air_properties_at_1_bar.dat", fu)
 26
 27   !-----------------------------!
 28   !   Read all the properties   !
 29   !-----------------------------!
 30   do i = 1, N_ITEMS
 31     call File % Read_Line(fu)
 32
 33     ! Read desired properties
 34     read(line % tokens(1), *) air_t(i)
 35     read(line % tokens(2), *) air_rho(i)
 36     read(line % tokens(3), *) air_mu(i)
 37     read(line % tokens(5), *) air_cp(i)
 38     read(line % tokens(6), *) air_lambda(i)
 39
 40     ! Fix units where needed (check the values in the table)
 41     air_cp(i) = air_cp(i) * 1.0e3
 42     air_mu(i) = air_mu(i) / 1.0e5
 43   end do
 44
 45   close(fu)
 46
 47   if(this_proc < 2) then
 48     print '(a)',        ' #============================================'
 49     print '(a)',        ' # Output from user function, read properties!'
 50     print '(a)',        ' #--------------------------------------------'
 51   end if
```

This probalby needs some explanation.  As arguments to this function, _Process_
sends a number of its classes which are used to model different aspects of numerical
simulation.  Class ```Field_Type``` holds velocities, temperatures and other variables
describing a flow _field_.  Class _Turb_Type_ holds variables describing the state
of turbulence; turbulent kinetic energy  (k), its dissipation (epsilon), individial
Reynolds stresses for second moment closures or turbulent statistics.  Multiphase
flows with interface tracking are described by class ```Vof_Type``` and Lagrangian
particle tracking with ```Swarm_Type```.  In addition to these classes, _Process_
also sends current time step (```n```) and physical time of simulation (```time```)

> **_Note:_** It may sound counter-intuitive to send the last two variables to a
procedure called at the beginning of simulation, but think of a restart.  This
procedure is also called, and time step and time of simulation are not zero any
more.

All the classes outlined above hold a pointer to a grid (```pnt_grid```) for which
are they defined. The grid and its entities could be accessed as ```Flow % pnt_grid```,
but we make sytnax shorter by introducing local pointer:
```
 15   type(Grid_Type), pointer :: Grid
```
and assingning it a value:
```
 20   Grid => Flow % pnt_grid
```

Line 20 calls T-Flows's class ```File_Type``` member function
```Open_For_Reading_Ascii``` whose purpose is clear from the name.  
The function returns a handle to file it openned called ```fu```.

From lines 30 to 43 we read the look-up table from the file, and store it into
memory defined in ```Types.f90```.  Note that in lines 41 and 42 we convert
units from the table to plain SI units for compatibility with T-Flows.  
While reading the files we use another procedure from ```File_Type``` called
```Read_Line``` which reads a line from ASCII file and splits it into individual
tokens.  Tokens are stored in the fields ```line % token(:)```.

Finally, in the lines 47 to 51 we print a message that physical properties
have been read.  Here we use global variable ```this_proc``` to make sure we
print the message only from first processor in parallel runs.

Once the look-up table is properly read into memory, we can use it at the
beginning of each time step with procedure ```Beginning_Of_Time_Step``` which reads:
```
  1 !==============================================================================!
  2   subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm, n, time)
  3 !------------------------------------------------------------------------------!
  4 !   This function is called at the beginning of time step.                     !
  5 !------------------------------------------------------------------------------!
  6   implicit none
  7 !---------------------------------[Arguments]----------------------------------!
  8   type(Field_Type),    target :: Flow
  9   type(Turb_Type),     target :: Turb
 10   type(Vof_Type),      target :: Vof
 11   type(Swarm_Type),    target :: Swarm
 12   integer, intent(in)         :: n     ! time step
 13   real,    intent(in)         :: time  ! physical time
 14 !-----------------------------------[Locals]-----------------------------------!
 15   type(Grid_Type), pointer :: Grid
 16   type(Var_Type),  pointer :: u, v, w, t, phi
 17   integer                  :: c, i
 18   real                     :: wi, wip
 19 !==============================================================================!
 20
 21   ! Take aliases
 22   Grid => Flow % pnt_grid
 23
 24   !------------------------------!
 25   !   Browse through all cells   !
 26   !------------------------------!
 27   do c = -Grid % n_bnd_cells, Grid % n_cells
 28
 29     ! Browse through all table entries
 30     do i = 1, N_ITEMS - 1
 31
 32       ! Did you find the right interval
 33       if(Flow % t % n(c) >= air_t(i) .and.  &
 34          Flow % t % n(c) <  air_t(i+1)) then
 35
 36         ! If so, calculate interpolation factors ...
 37         wi  = (air_t(i+1) - Flow % t % n(c)) / (air_t(i+1) - air_t(i))
 38         wip = 1.0 - wi
 39
 40         ! ... and interpolate physical properties
 41         Flow % density(c)      = wi * air_rho   (i)  + wip * air_rho   (i+1)
 42         Flow % viscosity(c)    = wi * air_mu    (i)  + wip * air_mu    (i+1)
 43         Flow % conductivity(c) = wi * air_lambda(i)  + wip * air_lambda(i+1)
 44         Flow % capacity(c)     = wi * air_cp    (i)  + wip * air_cp    (i+1)
 45       end if
 46
 47     end do
 48   end do
 49
 50   end subroutine
```
Arguments are the same as in the call to previous procedure and don't need explanation.  
What is new here are the fields ```n_bnd_cells``` and ```n_cells``` from class
```Grid_Type```..  They hold number of boundary cells and number of inside cells
respectivelly.  As you can see from the beginning of the loop in line 27,
boundary cells are stored with negative indices.  

Anyhow, for each of the cells in the grid, the entire look-up table is browsed
through in lines 30 - 47, and for temperature in that cell (```Flow % t % n(c)```)
the interval for interpolation is searched in lines 33 and 34.  Once found,
the interpolation factors are calculated in lines 37 and 38, and physical
properties interpolated in lines 41 - 44.

The only novelty in the present ```control``` file is that we set the time step
and the total number of time steps as:
```
 TIME_STEP                1.0e-2
 NUMBER_OF_TIME_STEPS  6000
```
resulting in a total physical time of simulation of 0.6 seconds.  

> **_Note:_** You might have noticed that we usuall set the number of time steps
and saving frequencies, to be multiples of 60.  We acknowledged that time steps
are usually set to be round in decimal values, such as 1.0e-3, 5e-4, ... and so
forth.  If one combines such values of time steps with number of time steps
which are multiples of 60, it is highly likely that simulation time and saving
frequency will easy to convert in minutes, maybe hours, which is easier to
comprehend than _thousands of seconds_.

There is nothing else particularly interseting in the ```control``` file for this case,
except the fact that physical properties are not defined.  It is because they are
set from user functions.

With all this explained, grid can be generated and converted with:
```
gmsh -3 air.geo
./Convert < convert.scr
./Process > out
```

> **_Note:_** Here we assume that you created soft links to executables with
script ```seek_binaries.sh``` and that you have compiled _Process_ with user
functions, i.e. with command: ```make DIR_CASE=../../Tests/Manual/Thermally_Driven/Varible```

In the case of thermally-driven cavity with Boussinesq approximation, we were
setting physical properties _and_ gravitational vector in a way to meet the
desired _Ra_ number.  Here, physical properties are physical, dependent on
temperature, as meassured and reported [here](https://theengineeringmindset.com/properties-of-air-at-atmospheric-pressure/).
So, one way to achive the desired _Ra_ would be to change boundary temperaures,
but that could put us out of the range of the look-up table we have.  The other
way, which was used here, was to change the charactersitic length and we do it
by setting a scaling factor in ```convert.scr```.  Please chech the ```readme```
file in the local directory for details and ```convert.scr``` in line 3 to see
the scaling factor used.

The final velocity vector imposed over the temperature fields look like this:

![Thermally-driven variable temp!](Documentation/Manual/Figures/thermally_driven_solution_variable_velo_temp.png "Thermally-driven variable temp")

whereas the same velocity vectors over pressure, look like:

![Thermally-driven variable press!](Documentation/Manual/Figures/thermally_driven_solution_variable_velo_press.png "Thermally-driven variable temp")

## Parallel processing  <a name="demo_parallel_proc"></a>

The cases considered in previous [chapter](#test_cases) were all on very small
grids an no parallel processing was needed.  As the grids we use grow bigger
and number of time steps we want to perform grow, the need for parallel runs
become imminent.

In this chapter, we will teach you:
- how to divide a computational grid over into a number of sub-domains
- how to conduct a parallel simulation

and, irrespective of parallel run or not:
- how to set up periodic cases with prescribed mass flux and pressure drops in
  periodic direction
- how to request from _Process_ to save results unplanned
- how to instruct _Process_ to exit prematurely.

Before you go any further, please make sure you have MPI installed on your system.  
The easiest way to check it is to try to run commands: ```mpirun``` and ```mpif90```.
If these commands exist, you are good to continue.  If not, you should install
MPI on your system.

> **_Note:_**  At the time of writing this manual, we did all developments
and testing with OpenMPI flavor of MPI, so it is a safer option.  Working with
MPICH shouldn't be an issue, either, but we are not sure at this moment.
Should you encounter problems with MPICH, feel free to contact the development
team.

### Compiling for parallel runs <a name="demo_parallel_proc_compiling"> </a>

If you are here, we believe that MPI is installed on your system and are
ready to continue.  Before starting a parallel run, you have to build the
executables.  The minimum you are going to need now, is:

1. _Generate_ or _Convert_,
2. _Divide_
3. _Process_

The _Generate_ or _Convert_ and _Divide_ are compiled in the usual way, that is
just by issuing command ```make``` in each of their sub-folder:
```
[root]/Sources/
├── Convert
├── Divide
├── Generate
└── Process
```

For the process, you still go to its folder but you should also tell ```make```
that you are building a parallel version of the code, hence command ```make MPI=yes```.
If you were compiling sequential version of the _Process_ before, you might
still have sequentially built objects in the local directory, so the safest
thing to do would be to issue these two commands in ```[root]/Sources/Process```:
```
make clean
make MPI=yes
```

> **_Warning:_** Please be reminded that you should compile them all of the
sub-programs with the same precision (```REAL=double```  or  ```REAL=single```)
option passed to ```make```.

### Creating and dividing the grid <a name="demo_parallel_proc_dividing"> </a>

We will use the case in ```[root]/Tests/Manual/Parallel``` to show you how to
set up a parallel run.  Please change to that folder and check the contents.  
They should read:
```
[root]/Tests/Manual/Parallel
├── control
├── convert.scr
└── rod_tet.geo
```

The purpose of these files should be clear to you by now.  If not, re-visit
[chapter](#test_cases).  The creation of the mesh is as usual:
```
gmsh -3 rod_tet.geo
```
followed by conversion process:
```
./Convert < convert.scr
```
After these steps, you should have the following files additional files in
the directory:
```
[root]/Tests/Manual/Parallel
├── ...
├── rod_tet.cfn
├── rod_tet.dim
├── rod_tet_dual.cfn
├── rod_tet_dual.dim
├── rod_tet_dual.edges.vtu
├── rod_tet_dual.faces.vtu
├── rod_tet_dual.shadows.vtu
├── rod_tet_dual.vtu
├── rod_tet.faces.vtu
├── ...
├── rod_tet.msh
├── ...
└── rod_tet.vtu
```

> **_Note:_** Only new files are shown here.  Old are represented with dots.

It should be clear to you by now that files ```.cfn``` and ```.dim``` is the grid
in internal T-Flows format.  The ```.vtu``` files come in two variants.  Without
insertion of ```_dual``` and with it, like we had in this [chapter](#demo_lid_driven_dual)
above.  That should tell you that the grid is polyhedral.  If you visualize
```rod_tet.vtu``` or ```rod_tet_dual.vtu```, the overall domain will look like
this:

![Rod bundle domain!](Documentation/Manual/Figures/rod_domain.png "Rod bundle domain")

It is a segment from a bundle of rods arranged in a hexagonal matrix.  We will
study a flow across this matrix with LES.  The details of the mesh for
```rod_tet.vtu``` or ```rod_tet_dual.vtu``` are different, as shown here for the former:

![Rod tet grid detail!](Documentation/Manual/Figures/rod_tet_detail.png "Rod tet grid detail")

and here for the latter:

![Rod tet dual grid detail!](Documentation/Manual/Figures/rod_tet_dual_detail.png "Rod tet dual grid detail")

Next step is to divide the mesh using the _Divide_.  The fastest is to invoke
_Divide_ from the command line as this:
```
./Divide  rod_tet_dual  6
```
by which we tell _Divide_ the name of the grid we want to divide, and number of
sub-divisions.  If the command is successful, your directory structure, showing
only new files, looks like this:
```
[root]/Tests/Manual/Parallel
├── ...
├── rod_tet_dual.pvtu
├── ...
├── Sub-00001
│   ├── rod_tet_dual.cfn
│   ├── rod_tet_dual.dim
│   └── rod_tet_dual.vtu
├── Sub-00002
│   ├── rod_tet_dual.cfn
│   ├── rod_tet_dual.dim
│   └── rod_tet_dual.vtu
├── Sub-00003
│   ├── rod_tet_dual.cfn
│   ├── rod_tet_dual.dim
│   └── rod_tet_dual.vtu
├── Sub-00004
│   ├── rod_tet_dual.cfn
│   ├── rod_tet_dual.dim
│   └── rod_tet_dual.vtu
├── Sub-00005
│   ├── rod_tet_dual.cfn
│   ├── rod_tet_dual.dim
│   └── rod_tet_dual.vtu
└── Sub-00006
    ├── rod_tet_dual.cfn
    ├── rod_tet_dual.dim
    └── rod_tet_dual.vtu
```
_Divide_ has created a number of sub-folder in the working directory called:
```Sub-00001``` to ```Sub-0000N```, where N is the number of processors (here
6 because that's what we requested.  Each of this sub-folder contains its
portion of the grid in T-Flows' format (the ```.cfn``` and ```.dim``` files)
and its portion of the ```.vtu``` file for visualization with ParaView.

Additional novelty is the file called ```rod_tet_dual.pvtu```, with extension
```.pvtu``` instead of ```.vtu```, where additional ```p``` is for _parallel_.
This file re-directs ParaView to read the partial ```.vtu``` files from the
new folder, and if you visualize the ```rod_tet_dual.pvtu``` it will show you
the domain decomposition obtained by _Divide_:

![Rod tet dual decomposed!](Documentation/Manual/Figures/rod_tet_dual_decomposed.png "Rod tet dual decomposed")

### Running the simulation in parallel <a name="demo_parallel_proc_running"> </a>

We are ready for parallel run.  Provided that you compiled the _Process_ with
```MPI=yes``` option, you can start the parallel run with ```mpirun``` command
as follows:
```
mpirun  -np 6  ./Process > out_6_proc  &
```
where ```-np 6``` tells ```mpirun``` how many processors you want to use.  And
that is it, the code should be running now.

#### Setting the mass fluxes and pressure drops

While the process is running, we should spare a few words on a couple of new
concepts introduced in the ```control``` file.  The flow across these rod
bundle is driven by a pressure drop in _x_ direction, which is re-computed
by T-Flows at every time step to give the desired mass flux in _x_ direction.
This is achieved with following two lines:
```
 PRESSURE_DROPS         3.6  0.0  0.0
 MASS_FLOW_RATES        0.5  0.0  0.0
```
If only first (```PRESSURE_DROPS```) was applied, T-Flows would keep it constant.
With second line (```MASS_FLOW_RATES```) it re-adjusts the pressure drop to
keep mass flux at the desired rate.  

> **_Note:_** The first option only tends to lead to a rather slow flow
development.

The mass flow rates and the pressure drops are printed by _Process_ at the end
of each time step, see the excerpt from the log file ```out_6_proc``` with
focus on the fields ```Flux x``` and ```Pdrop x```:
```
  #=============================================================================================================================#
  #                                                        Iteration:  2                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  1 4.883E-04 | V   :  1 3.292E-04 | W   :  1 2.223E-06 | PP  : 58 8.465E-07 | MASS:    3.230E-08 |                    #
  #=============================================================================================================================#
  #                                                        Iteration:  3                                                        #
  #--------------------+--------------------+--------------------+--------------------+--------------------+--------------------#
  # U   :  1 2.440E-04 | V   :  1 1.513E-04 | W   :  1 1.432E-06 | PP  : 58 8.293E-07 | MASS:    2.242E-08 |                    #
  #--------------------+=========================================+=========================================+--------------------#
                       #    Maximum Courant number: 6.114E-01    |    Maximum Peclet number: 2.062E+03     #
                       #---------------------------+-------------+-------------+---------------------------#
                       #    Flux x :  5.000E-01    |    Flux y :  -1.15E-05    |    Flux z :  -2.19E-09    #
                       #    Pdrop x:  3.601E-02    |    Pdrop y:   0.00E+00    |    Pdrop z:   0.00E+00    #
                       #---------------------------+---------------------------+---------------------------#
```
The mass fluxes in question are computed in three characteristic planes (orthogonal
to _x_, _y_ and _z_ axis at the point defined in the control file as:
```
POINT_FOR_MONITORING_PLANES    0.007  0.00  0.002
```
For flows with prescribed mass fluxes in periodic direction such as this one,
it is a good practice to define this point.

#### Saving and/or exiting prematurely.

For this case, we set the desired number of time steps to 6000, and we set
ssaving interval to each 1200 time steps.  We know it is a big grid, and we
don't want to overfill the disk.  Now imagine that you curious to see the
results before time step reaches the prescribed interval.  (That is, in
essence, not a bad idea, as you can use to make sure results are not marred
with numerical instabilities.)

In order to achieve that, you should create a file called ```save_now``` in
the directory where the simulation is running.  The easiest way to do it on
Linux is to issue the ```touch``` command like this:
```
touch save_now
```
If _Process_ finds this file in the working directory, it will save the results
immediately, delete the file ```save_now``` and continue the simulation.

If we do it now, at the time of writing this manual, _Process_ will save results
at the time step of ... 357.  File ```rod_tet_dual-ts000357.pvtu``` has been
created and we can visualize it to show current velocity field in _x_ direction:

![Rod tet dual vel x!](Documentation/Manual/Figures/rod_tet_dual_ts0357_vel_x.png "Rod tet dual")

This actually happened to be a lucky moment, because we can observe recirculation
zones behind each of the rods (blue areas).  A more detailed view of the recirculation
zones shown as vector imposed over pressure looks like:

![Rod tet dual vel x!](Documentation/Manual/Figures/rod_tet_dual_ts0357_recirc.png "Rod tet dual")

Next time we _touched_ ```save_now``` happened to be time step 1061, at which the
flow start to exhibit three-dimensional patterns:

![Rod tet dual vel x!](Documentation/Manual/Figures/rod_tet_dual_ts1061_vel_x.png "Rod tet dual")

Fine, maybe the simulation doesn't need to continue.  You have hopefully grasped
how to launch a parallel simulation.  Before we end, let's just stop the simulation
which is running.  We could do it with a command:
```
pkill Process
```
which is a rather abrupt way to stop, but we may also direct _Process_ to end
gracefully, saving the last results and backup files before it stops.  To do
that, we create a file called ```exit_now``` in the running directory.  If
_Process_ finds this file, it will save results and backup, and exit.  Before
exiting, it will delete the ```exit_now``` file to prevent interfering with the
next run in the same directory.

#### Cleaning the directory from ```.vtu``` files

Parallel runs, almost invariantly, create a huge amount of data on the file
system.  If you want to clear a certain folder in which a parallel simulation
was running in order to free some disk space, deleting ```.pvtu``` files
alone will be of much benefit.  The sub-folders holding sub-domains
(```Sub-00001``` - ```Sub-0000N```) are the biggest in volume.  So, in order
to get rid of ```.vtu``` files, you should delete like this:
```
rm -f Sub-0*/*.vtu
```

#### Thing to try next

Another interesting feature of this case is that it has no less than four
periodic direction.  Instead of converting the grid obtaing from GMSH with:
```
./Convert < convert.scr
```
try to do it without the script, simply by:
```
./Convert
```
For this case, it is interestng to see how the periodic directions are
dealt with.

# Benchmark cases  <a name="bench_cases"></a>

## Large eddy simulation over a matrix of cubes] <a name="bench_cases_matrix"> </a>
