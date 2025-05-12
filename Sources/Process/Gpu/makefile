#-------------------------------------------------------------------------------
#
#   Process Makefile
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

#-------------------------------------
#   Check values passed to makefile
#-------------------------------------
include makefile_check_arguments

# MPI in addition
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

# GPU in addition
ifneq ($(GPU), yes)
  ifneq ($(GPU), no)
    $(info #=================================================================)
    $(info # Error:                                                          )
    $(info #   Possible values for GPU are "yes" or "no" and they are        )
    $(info #   case sensitive.  You used the "$(GPU)", which is unknown.     )
    $(info #-----------------------------------------------------------------)
    $(error This error is critical, exiting!)
  endif
endif

#--------------------------
#
#   Variable definitions
#
#--------------------------

#---------------------------------------------------------------
#   Directories for objects and modules. (No need to change.)
#---------------------------------------------------------------
DIR_BINARY = ../../../Binaries/Gpu
DIR_SHARED = ../../Shared
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
DIR_USER   = User_Mod
DIR_CASE   = no
VPATH = $(DIR_SHARED):.

#---------------------------------------------------
#   Dissengage GPU for Intel and Gnu and for OpenMP
#---------------------------------------------------
ifeq ($(FORTRAN), intel)
  GPU = no
endif

ifeq ($(FORTRAN), gnu)
  GPU = no
endif

ifeq ($(OMP), yes)
  GPU = no
endif

#------------------------------
#   Python for preprocessing
#------------------------------
PYTHON       = python3
PREPROCESSOR = tf-acc.py
PREOPTION    = openacc

ifeq ($(OMP), yes)
  PREOPTION = openmp
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
$(info #        <DEBUG=no/yes> <ASSERT=yes/no> <SORT=quick/heap> <MPI=no/yes>  )
$(info #        <SORT=quick/heap> <DIR_CASE=full_or_relative_path_to_case>     )
$(info #        <clean/preprocess>                                             )
$(info #                                                                       )
$(info # Notes:                                                                )
$(info #   1. The first item, for each of the options above, is the default.   )
$(info #   2. If you use gnu or intel compiler, GPU will be set to "no"        )
$(info #   3. Target "clean" removes all compiled files and binaries.          )
$(info #   4. Target "preprocess" will only preprocess .fpp to .f90 files.     )
$(info #      Effectivelly, this targets inserts OpenACC into Fortran sources. )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make              - compile with nvidia compiler with GPU           )
$(info #   make GPU=no       - compile with nvidia compiler without GPU        )
$(info #   make FORTAN=intel - compile with intel compiler                     )
$(info #   make DEBUG=yes    - compile with nvidia compiler in debug mode      )
$(info #   make clean        - remove all compiled files and binaries          )
$(info #-----------------------------------------------------------------------)

#-------------------
#   METIS library
#-------------------
METIS_LIB = ../../Libraries/Metis_5.1.0_Linux_64/libmetis_i32_r64.a
ifeq ($(REAL), single)
  METIS_LIB = ../../Libraries/Metis_5.1.0_Linux_64/libmetis_i32_r32.a
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
    OPT_F_COMP += -O2
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

  # Disable suggestion to use new version of Fortran compiler
  OPT_F_COMP += -diag-disable=10448

  $(info # Using Intel Fortran compiler with options: $(OPT_F_COMP))

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
  ifeq ($(GPU), yes)
    OPT_F_COMP += -acc -Minfo=accel
  endif
  ifeq ($(OMP), yes)
    OPT_F_COMP += -mp
  endif

  $(info # Using Nvidia Fortran compiler with options: $(OPT_F_COMP))

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
   OPT_F_COMP += -O0 -g -fcheck=all -Wall -Wno-integer-division
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

  $(info # Using GNU Fortran compiler with options: $(OPT_F_COMP))
endif

# Set initial options for C compiler - this is probably not needed
OPT_C_COMP = -cpp

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

#--------------------
#   User functions
#--------------------
ifneq ($(DIR_CASE), no)
  ifeq ('$(wildcard $(DIR_CASE))','')
    $(info Case directory $(DIR_CASE) does not exist!)
  else
    $(info Compling for the case in $(DIR_CASE))
    ifeq ('$(wildcard $(DIR_CASE)/$(DIR_USER))','')
      $(info Directory with user functions $(DIR_CASE)/$(DIR_USER) does not exist!)
    else
      $(info Linking user functions from directory $(DIR_CASE)/$(DIR_USER))

      # Relative path to $DIR_CASE from ./
      REL_DIR=$(shell realpath --relative-to="./" "$(DIR_CASE)")

      # Fetch all user functions from the case directory
      STEP_A  = $(wildcard $(REL_DIR)/$(DIR_USER)/*.f90) $(wildcard $(REL_DIR)/$(DIR_USER)/*.fpp)
      STEP_B  = $(notdir $(STEP_A))

      # Append "../" to each of the file names
      SRC_USR = $(addprefix ../,$(STEP_A))
      $(info $(SRC_USR))

      ifneq ("$(wildcard $(STEP_A))","")
        # Make soft links
        $(shell ln -f -s $(SRC_USR) $(DIR_USER))
      endif

    endif
  endif
endif

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
             Gpu_Mod.f90            \
             Turb_Mod.f90           \
             User_Mod.f90           \
             Backup_Mod.f90         \
             Results_Mod.f90        \
             Read_Controls_Mod.f90  \
             Process_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all functions are obtained by a shell command
SRC_F_FUN = $(shell ls -1 *.f90                        \
                        | xargs -n1 basename   \
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
#   Python Preprocessing Rule for .fpp files
#---------------------------------------------------------
# This rule preprocesses .fpp files into .f90 files
%.f90: %.fpp
	@echo "Preprocessing $< with Python to generate $@"
	@$(PYTHON) $(PREPROCESSOR) $(PREOPTION) $< $@

# Add the .fpp files here, which need preprocessing to .f90
SRC_FPP = $(shell find . -name "*.fpp")

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
	@echo -n Removing objects and modules ...
	find . -name "*.fpp" | xargs touch
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
	@echo done
	@echo -n Removing links to sources in case directory ...
	@rm -f $(DIR_USER)/*.f??
	@echo done
	@echo -n Retreiving default user functions from git repository ...
	@git checkout $(DIR_USER)/*.f??
	@echo done

#---------------------------------------------------------
#   Preprocessing only target
#---------------------------------------------------------
# This target will only preprocess the .fpp files
preprocess: $(SRC_FPP:.fpp=.f90)
	@echo "Preprocessing complete."

