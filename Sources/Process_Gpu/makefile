#===============================================================================
#
#   SFS Makefile
#
#-------------------------------------------------------------------------------

#-------------------------------------
#   Default options for compilation
#-------------------------------------
FORTRAN ?= nvidia
GPU     ?= yes
ASSERT  ?= yes
DEBUG   ?= no
PROF    ?= no
REAL    ?= double
MPI     ?= no
OMP     ?= no
SORT    ?= quick

#--------------------------
#   Variable definitions
#--------------------------

# MPI
ifneq ($(MPI), yes)
  ifneq ($(MPI), no)
    $(info #=================================================================)
    $(info # Error:                                                          )
    $(info #   Possible values for MPI are "yes" or "no" and they are        )
    $(info #   case sensitive.  You used the "$(MPI)", which is unknown.     )
    $(info #-----------------------------------------------------------------)
    $(error This error is critical, exiting!)
  endif
endif

#---------------------------------------------------------------
#   Directories for objects and modules. (No need to change.)
#---------------------------------------------------------------
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
VPATH = $(DIR_SHARED):.

#--------------------------------------
#   Dissengage GPU for Intel and Gnu
#--------------------------------------
ifeq ($(FORTRAN), intel)
  GPU = no
endif

ifeq ($(FORTRAN), gnu)
  GPU = no
endif

# Program name (This should hardly change)
PROGRAM_NAME = Process_Cpu
ifeq ($(GPU), yes)
  PROGRAM_NAME = Process_Gpu
endif

PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

#----------------------------------------------------------
#   Initialize variables to pass on to preprocessor with
#   variable T_FLOWS_COMPILATION which, if set to one,
#   means unit is compiled from T-Flow
#   It later adds T_FLOWS_MPI and T_FLOWS_PETSC
#----------------------------------------------------------
PASS_ON = -DT_FLOWS_COMPILATION=1

#----------------------------------------------------------
#   Add T_FLOWS_PROGRAM, T_FLOWS_ASSERT and T_FLOWS_DEBUG
#   and later on it add also T_FLOWS_MPI and T_FLOWS_PETSC
#----------------------------------------------------------
PASS_ON += -DT_FLOWS_PROGRAM=4

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

ifeq ($(GPU), yes)
  PASS_ON += -DT_FLOWS_GPU=1
else
  PASS_ON += -DT_FLOWS_GPU=0
endif

$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with $(FORTRAN) compiler                    )
$(info #-----------------------------------------------------------------------)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=nvidia/gnu/intel> <GPU=yes/no> <REAL=double/single>   )
$(info #        <DEBUG=no/yes> <ASSERT=yes/no> <SORT=quick/heap> <OMP=no/yes>  )
$(info #                                                                       )
$(info # Notes:                                                                )
$(info #   1. The first item, for each of the options above, is the default.   )
$(info #   2. If you use gnu or intel compiler, GPU will be set to "no"        )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make              - compile with nvidia compiler with GPU           )
$(info #   make GPU=no       - compile with nvidia compiler without GPU        )
$(info #   make FORTAN=intel - compile with intel compiler                     )
$(info #   make DEBUG=yes    - compile with nvidia compiler in debug mode      )
$(info #-----------------------------------------------------------------------)

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
  FC = ifort
  ifeq ($(DEBUG),yes)
    OPT_F_COMP = -module $(DIR_MODULE) -cpp -r8 -O0 -g -warn all -check all \
               -debug all -fpe-all=0 -traceback
  else
    OPT_F_COMP = -module $(DIR_MODULE) -cpp -r8 -O3
  endif

  # Disable suggestion to use new version of Fortran compiler
  OPT_F_COMP += -diag-disable=10448

  $(info # Using Intel Fortran compiler with options: $(OPT_F_COMP))
endif

#------------
#   Nvidia
#------------

# Fortran == nvidia
ifeq ($(FORTRAN), nvidia)
  FC = nvfortran
  ifeq ($(DEBUG),yes)
    OPT_F_COMP = -module $(DIR_MODULE) -cpp -r8 -O0 -g
  else
    OPT_F_COMP = -module $(DIR_MODULE) -cpp -r8 -O3
  endif
  ifeq ($(GPU),yes)
    OPT_F_COMP += -acc -Minfo=accel
  endif

  $(info # Using Nvidia Fortran compiler with options: $(OPT_F_COMP))
endif

#---------
#   Gnu
#---------

# Fortran == gnu
ifeq ($(FORTRAN), gnu)
  FC = gfortran
  ifeq ($(DEBUG),yes)
    OPT_F_COMP = -J $(DIR_MODULE) -cpp -fdefault-real-8 -O0 -g  \
                                  -ffree-line-length-none       \
                                  -Wunused-parameter -Warray-temporaries -Wall
  else
    OPT_F_COMP = -J $(DIR_MODULE) -cpp -fdefault-real-8 -O3     \
                                  -ffree-line-length-none
  endif

  $(info # Using GNU Fortran compiler with options: $(OPT_F_COMP))
endif

#-----------------------------------------------------------
#   Set linker options to be the same as compiler options
#-----------------------------------------------------------
OPT_LINK = $(OPT_F_COMP)

#-----------------------
#   MPI related stuff
#-----------------------

# Use here your command for mpi fortran compiler:
# (This is actually not too good, because not all
# compiler vendors have MPI wrapper named mpif90)
ifeq ($(MPI), yes)
  ifeq ($(FCOMP),)
    ifeq ($(FORTRAN), intel)
      FC = mpiifx
    else
      FC = mpif90
    endif
  else
    FC = $(FCOMP)
  endif
endif

# Variable T_FLOWS_MPI is passed on to compilers
ifeq ($(MPI), yes)
  PASS_ON += -DT_FLOWS_MPI=1
else
  PASS_ON += -DT_FLOWS_MPI=0
endif

#------------------------------------------------------
#   List of sources for modules and functions
#------------------------------------------------------
#   Modules' order must obey their dependency
#   This list should therefore be written "by hand".
#------------------------------------------------------

#-------------
#   Modules
#-------------

# Modules in shared directory
SRC_F_MOD = Const_Mod.f90         \
            Region_Mod.f90        \
            String_Mod.f90        \
            Comm_Mod.f90          \
            Assert_Mod.f90        \
            Tokenizer_Mod.f90     \
            Message_Mod.f90       \
            Memory_Mod.f90        \
            Omp_Mod.f90           \
            Math_Mod.f90          \
            Swap_Mod.f90          \
            Sort_Mod.f90          \
            File_Mod.f90          \
            Control_Mod.f90       \
            Profiler_Mod.f90      \
            Vtk_Mod.f90           \
            Metis_Mod.f90         \
            Grid_Mod.f90          \
            Work_Mod.f90

# Modules in local directory
SRC_F_MOD += Gpu_Pointers_Mod.f90   \
             Bulk_Mod.f90           \
             Face_Mod.f90           \
             Sparse_Mod.f90         \
             Linalg_Mod.f90         \
             Native_Mod.f90         \
             Numerics_Mod.f90       \
             Var_Mod.f90            \
             Time_Mod.f90           \
             Field_Mod.f90          \
             Iter_Mod.f90           \
             Info_Mod.f90           \
             Turb_Mod.f90           \
             Gpu_Mod.f90            \
             Backup_Mod.f90         \
             Results_Mod.f90        \
             Read_Controls_Mod.f90  \
             Process_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all functions are obtained by a shell command
SRC_F_FUN = $(shell ls -1 *.f90			\
                        | xargs -n1 basename	\
                        | grep -v -i _Mod)

# Sources for all shared functions
SRC_F_FUN += Key_Ind.f90

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

# Modules
$(DIR_OBJECT)/%.o: %.f90 %/*.f90 makefile*
	@echo $(FC) $<
	@$(FC) $(OPT_F_COMP) $(PASS_ON) -c -o $@ $<

# Functions
$(DIR_OBJECT)/%.o: %.f90 makefile*
	@echo $(FC) $<
	@$(FC) $(OPT_F_COMP) $(PASS_ON) -c -o $@ $<

#-----------------------------------
#   Rule to build main program
#-----------------------------------
#   Note: Should not be modified.
#-----------------------------------
$(PROGRAM_FILE): $(OBJ)
	@echo Linking "\033[0;32m $(METIS_LIB) $(PROGRAM_FILE) \033[0m"
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
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
