#===============================================================================
#
#   Generate Makefile
#
#-------------------------------------------------------------------------------

#-------------------------------------
#   Default options for compilation
#-------------------------------------
FORTRAN ?= gnu
ASSERT  ?= yes
DEBUG   ?= no
PROF    ?= no
REAL    ?= double
OMP     ?= yes
SORT    ?= quick

#-------------------------------------
#   Check values passed to makefile
#-------------------------------------
include makefile_check_arguments

#---------------------------------------------------------------
#   Directories for objects and modules. (No need to change.)
#---------------------------------------------------------------
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
VPATH = $(DIR_SHARED):.

# Program name (This should hardly change)
PROGRAM_NAME = Generate
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

#----------------------------------------------------------
#   Initialize variables to pass on to preprocessor with
#   variable T_FLOWS_COMPILATION which, if set to one,
#   means unit is compiled from T-Flow
#----------------------------------------------------------
PASS_ON = -DT_FLOWS_COMPILATION=1

#----------------------------------------------------------
#   Add T_FLOWS_PROGRAM, T_FLOWS_ASSERT and T_FLOWS_DEBUG
#   and later on it add also T_FLOWS_MPI and T_FLOWS_PETSC
#----------------------------------------------------------
PASS_ON += -DT_FLOWS_PROGRAM=2

ifeq ($(ASSERT), yes)
  PASS_ON += -DT_FLOWS_ASSERT=1
else
  PASS_ON += -DT_FLOWS_ASSERT=0
endif

ifeq ($(DEBUG), yes)
  PASS_ON += -DT_FLOWS_DEBUG=1
else
  PASS_ON += -DT_FLOWS_DEBUG=0
endif

ifeq ($(SORT), quick)
  PASS_ON += -DT_FLOWS_QUICKSORT=1
else
  PASS_ON += -DT_FLOWS_QUICKSORT=0
endif

$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with $(FORTRAN) compiler                    )
$(info #-----------------------------------------------------------------------)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=gnu/intel/nvidia> <DEBUG=no/yes> <ASSERT=yes/no>      )
$(info #        <FCOMP=gfortran/ifort/nvfortran/mpif90/mpiifort/...>           )
$(info #        <PROF=no/yes> <REAL=double/single> <OMP=yes/no>                )
$(info #        <SORT=quick/heap> <clean>                                      )
$(info #                                                                       )
$(info # Notes:                                                                )
$(info #   1. The first item, for each of the options above, is the default.   )
$(info #   2. Target "clean" removes all compiled files and binaries.          )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make               - compile with gnu compiler                      )
$(info #   make FORTRAN=intel - compile with intel compiler                    )
$(info #   make DEBUG=yes     - compile with gnu compiler in debug mode        )
$(info #   make ASSERT=no     - compile without assert with gnu compiler       )
$(info #   make OMP=yes       - compile with gnu compiler for parallel Generate)
$(info #   make SORT=heap     - use heap-sort variant of sorting functions     )
$(info #   make clean         - remove all compiled files and binaries         )
$(info #-----------------------------------------------------------------------)

#----------------------------------------------------
# Error trap for unknow options passed to makefile
#----------------------------------------------------
KNOWN_VARIABLES := FORTRAN FCOMP ASSERT DEBUG PROF REAL OMP SORT
COMMAND_LINE_VARIABLES = $(foreach v,$(.VARIABLES), \
  $(if $(findstring command line,$(origin $(v))),$(v),))
UNKNOWN_VARIABLES = $(filter-out $(KNOWN_VARIABLES),$(COMMAND_LINE_VARIABLES))
ifneq ($(strip $(UNKNOWN_VARIABLES)),)
$(info # Error:                                                                )
$(info #   Some options you passed are misspelled or not known to makefile.    )
$(info #   Check the list of valid options above, in the section "Usage".      )
$(info #   Unknown options: $(UNKNOWN_VARIABLES)                               )
$(info #-----------------------------------------------------------------------)
$(error This error is critical, exiting!)
endif

#-------------------
#   METIS library
#-------------------
METIS_LIB = ../Libraries/Metis_5.1.0_Linux_64/libmetis_i32_r64.a
ifeq ($(REAL), single)
  METIS_LIB = ../Libraries/Metis_5.1.0_Linux_64/libmetis_i32_r32.a
endif

#-------------------------------------------------------------------------
#
#   Compiler and linker options
#
#-------------------------------------------------------------------------
#   Note: Changes only when support to a new Fortran compiler is added.
#-------------------------------------------------------------------------

#-----------
#   Intel
#-----------

# Fortran == intel
ifeq ($(FORTRAN), intel)
  ifeq ($(FCOMP),)
    FC = ifort
  else
    FC = $(FCOMP)
  endif
  OPT_F_COMP = -module $(DIR_MODULE) -cpp
  ifeq ($(DEBUG), yes)
    OPT_F_COMP += -O0 -g -traceback -check bounds
  else
    OPT_F_COMP += -O3
  endif
  ifeq ($(PROF), yes)
    OPT_F_COMP += -g -O2 -shared-intel -debug inline-debug-info  \
                  -D TBB_USE_THREADING_TOOLS
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -r8
  endif
  ifeq ($(OMP), yes)
    OPT_F_COMP += -qopenmp
    ifeq ($(PROF), yes)
      OPT_F_COMP += -qopenmp-link dynamic -parallel-source-info=2
    endif
  endif

#------------
#   Nvidia
#------------

# Fortran == nvidia
else ifeq ($(FORTRAN), nvidia)
  ifeq ($(FCOMP),)
    FC = nvfortran
  else
    FC = $(FCOMP)
  endif
  OPT_F_COMP = -module $(DIR_MODULE) -cpp -Mbackslash
  ifeq ($(DEBUG), yes)
    OPT_F_COMP += -O0 -g
  else
    OPT_F_COMP += -O3
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -r8
  endif
  ifeq ($(OMP), yes)
    OPT_F_COMP += -mp
  endif

#---------
#   Gnu
#---------

# Fortran == gnu
else
  ifeq ($(FCOMP),)
    FC = gfortran
  else
    FC = $(FCOMP)
  endif
  OPT_F_COMP = -J $(DIR_MODULE) -cpp -ffree-line-length-none
  ifeq ($(DEBUG), yes)
    OPT_F_COMP += -O0 -g -fcheck=all -Wunused-parameter -Wall  \
                  -Wno-integer-division
  else
    OPT_F_COMP += -O3
  endif
  ifeq ($(PROF), yes)
    OPT_F_COMP += -pg -g
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -fdefault-real-8
  endif
  ifeq ($(OMP), yes)
    OPT_F_COMP += -fopenmp
  endif
endif

#-----------------------------------------------------------
#   Set linker options to be the same as compiler options
#-----------------------------------------------------------
OPT_LINK = $(OPT_F_COMP)

#----------------------------------------------------------------
#   MPI related stuff (Generate is not an MPI program, pass 0)
#----------------------------------------------------------------
PASS_ON += -DT_FLOWS_MPI=0

#------------------------------------------------------
#   List of sources for modules and functions
#------------------------------------------------------
#   Modules' order must obey their dependency
#   This list should therefore be written "by hand".
#------------------------------------------------------

#-------------
#   Modules
#-------------

# Modules in shared directories
SRC_F_MOD = Const_Mod.f90      \
            Region_Mod.f90     \
            String_Mod.f90     \
            Comm_Mod.f90       \
            Assert_Mod.f90     \
            Tokenizer_Mod.f90  \
            Message_Mod.f90    \
            Memory_Mod.f90     \
            Control_Mod.f90    \
            Profiler_Mod.f90   \
            Math_Mod.f90       \
            File_Mod.f90       \
            Vtk_Mod.f90        \
            Metis_Mod.f90      \
            Swap_Mod.f90       \
            Sort_Mod.f90       \
            Grid_Mod.f90

# Modules in the local directory
SRC_F_MOD += Point_Mod.f90     \
             Line_Mod.f90      \
             Block_Mod.f90     \
             Range_Mod.f90     \
             Smooths_Mod.f90   \
             Refines_Mod.f90   \
             Domain_Mod.f90    \
             Generate_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all shared functions
SRC_F_FUN = Probe_1d_Nodes.f90

# Local functions are added one by one ...
SRC_F_FUN += Main_Gen.f90

#----------------------------------------------------------------------
#   List of objects generated from the list of modules and functions
#----------------------------------------------------------------------
#   Note: This doesn't need editing.
#----------------------------------------------------------------------
OBJ_F_MOD = $(SRC_F_MOD:%.f90=$(DIR_OBJECT)/%.o)
OBJ_F_FUN = $(SRC_F_FUN:%.f90=$(DIR_OBJECT)/%.o)
OBJ = $(OBJ_F_MOD) $(OBJ_F_FUN)

#---------------------------------------------------------
#   Default rule to build Fortran modules and functions
#---------------------------------------------------------
#   Note: This doesn't need editing.
#---------------------------------------------------------

# Fortran sources
$(DIR_OBJECT)/%.o: %.f90
	@echo $(FC) $<
	@$(FC) $(OPT_F_COMP) $(PASS_ON) -c -o $@ $<

#-----------------------------------
#   Rule to build main program
#-----------------------------------
#   Note: Should not be modified.
#-----------------------------------
$(PROGRAM_FILE): $(OBJ)
	@echo -e Linking "\033[0;32m $(METIS_LIB) $(PROGRAM_FILE) \033[0m"
	@$(FC) $(OPT_LINK) -o $(PROGRAM_FILE) $(OBJ) $(METIS_LIB)

#--------------------------------------------------------------------
#   Explicit dependencies for modules
#--------------------------------------------------------------------
#   These are automatically generated by:
#   Sources/Utilities/create_external_dependencies_for_makefile.sh
#--------------------------------------------------------------------
include makefile_explicit_dependencies

#---------------------
#   Explicit target
#---------------------
clean:
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod \
	      $(DIR_MODULE)/*.f90 $(PROGRAM_FILE)
