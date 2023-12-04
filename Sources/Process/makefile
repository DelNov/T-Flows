#-------------------------------------------------------------------------------
#
#   Process Makefile
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
MPI     ?= no
OMP     ?= no
SORT    ?= quick

#-------------------------------------
#   Check values passed to makefile
#-------------------------------------
include makefile_check_arguments

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
DIR_USER   = User_Mod
DIR_CASE   = no
VPATH = $(DIR_SHARED):.

# Program name (This should hardly change)
PROGRAM_NAME = Process
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

$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with $(FORTRAN) compiler                    )
$(info #-----------------------------------------------------------------------)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=gnu/intel/nvidia> <DEBUG=no/yes> <ASSERT=yes/no>      )
$(info #        <FCOMP=gfortran/ifort/nvfortran/mpif90/mpiifort/...>           )
$(info #        <PROF=no/yes> <REAL=double/single> <MPI=no/yes> <OMP=no/yes>   )
$(info #        <SORT=quick> <DIR_CASE=full_or_relative_path_to_case>          )
$(info #                                                                       )
$(info # Notes:                                                                )
$(info #   1. The first item, for each of the options above, is the default.   )
$(info #   2. If PETSC_DIR environment variable is defined, the code will      )
$(info #      make an attempt to compile with PETSc libraries.                 )
$(info #   3. If PETSC_DIR and PETSC_ARCH set the environment with different   )
$(info #      compiler than you are trying to use, the compilation might fail. )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make                        - sequential version with gnu compiler  )
$(info #   make MPI=yes                - MPI version with gnu compiler         )
$(info #   make MPI=yes FCOMP=compiler - MPI version with custom compiler      )
$(info #   make FORTRAN=intel MPI=yes  - MPI version with intel compiler       )
$(info #   make MPI=yes SORT=heap      - MPI version with gnu compiler and     )
$(info #                               - heap-sort version of sorting          )
ifneq ($(PETSC_DIR),)
$(info #                                                                       )
$(info # Warning:                                                              )
$(info #   At the very end of the makefile, environment variables for PETSc    )
$(info #   libraries are included.  If you have a compiler mismatch there,     )
$(info #   that is if you PETSc was compiled with different compiler from      )
$(info #   the one you are using to compile Process, errors will occur.        )
endif
$(info #-----------------------------------------------------------------------)

#----------------------------------------------------
# Error trap for unknow options passed to makefile
#----------------------------------------------------
KNOWN_VARIABLES := FORTRAN FCOMP ASSERT DEBUG PROF REAL MPI OMP DIR_CASE SORT
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
  OPT_F_COMP_USER = $(OPT_F_COMP)

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
  OPT_F_COMP_USER = $(OPT_F_COMP)

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
  OPT_F_COMP_USER = $(OPT_F_COMP)
  OPT_F_COMP_USER += -Wno-unused-dummy-argument
  OPT_F_COMP_USER += -Wno-unused-variable

endif

# Set initial options for C compiler - this is probably not needed
OPT_C_COMP = -cpp

#----------------------------------------------
#   Additional options to compile with PETSc
#----------------------------------------------
ifneq ($(PETSC_DIR),)
  PASS_ON += -DT_FLOWS_PETSC=1
  OPT_F_COMP += -I ${PETSC_DIR}/include
  OPT_F_COMP += -I ${PETSC_DIR}/${PETSC_ARCH}/include
  OPT_C_COMP += -I ${PETSC_DIR}/include
  OPT_C_COMP += -I ${PETSC_DIR}/${PETSC_ARCH}/include
else
  PASS_ON += -DT_FLOWS_PETSC=0
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
      FC = mpiifort
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
      STEP_A  = $(shell ls -1 $(REL_DIR)/$(DIR_USER)/*.f90)
      STEP_B  = $(notdir $(STEP_A))

      # Append "../" to each of the file names
      SRC_USR = $(patsubst %.f90, ../%.f90, $(STEP_A))
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
            Memory_Mod.f90        \
            Region_Mod.f90        \
            String_Mod.f90        \
            Comm_Mod.f90          \
            Assert_Mod.f90        \
            Tokenizer_Mod.f90     \
            Message_Mod.f90       \
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
            Work_Mod.f90          \
            Stl_Mod.f90           \
            Iso_Polygons_Mod.f90  \
            Polyhedron_Mod.f90    \
            Isoap_Mod.f90

# Sources for modules in the local directory
SRC_F_MOD += Point_Mod.f90          \
             Bulk_Mod.f90           \
             Face_Mod.f90           \
             Vert_Mod.f90           \
             Elem_Mod.f90           \
             Side_Mod.f90           \
             Matrix_Mod.f90         \
             Vector_Mod.f90         \
             Native_Mod.f90         \
             Petsc_Mod.f90          \
             Var_Mod.f90            \
             Solver_Mod.f90         \
             Time_Mod.f90           \
             Numerics_Mod.f90       \
             Field_Mod.f90          \
             Iter_Mod.f90           \
             Front_Mod.f90          \
             Surf_Mod.f90           \
             Info_Mod.f90           \
             Eddies_Mod.f90         \
             Monitor_Mod.f90        \
             Porosity_Mod.f90       \
             Turb_Mod.f90           \
             Vof_Mod.f90            \
             Particle_Mod.f90       \
             Swarm_Mod.f90          \
             Interface_Mod.f90      \
             User_Mod.f90           \
             Backup_Mod.f90         \
             Results_Mod.f90        \
             Read_Controls_Mod.f90  \
             Process_Mod.f90

SRC_C_FUN = C_Petsc_Interface.c

#---------------
#   Functions
#---------------

# Sources for all functions in the local directory
SRC_F_FUN = Main_Pro.f90

# Sources for all shared functions
SRC_F_FUN += Key_Ind.f90

#----------------------------------------------------------------------
#   List of objects generated from the list of modules and functions
#----------------------------------------------------------------------
#   Note: This doesn't need editing.
#----------------------------------------------------------------------
OBJ_F_MOD = $(SRC_F_MOD:%.f90=$(DIR_OBJECT)/%.o)
OBJ_F_FUN = $(SRC_F_FUN:%.f90=$(DIR_OBJECT)/%.o)
OBJ_C_FUN = $(SRC_C_FUN:%.c=$(DIR_OBJECT)/%.o)
OBJ = $(OBJ_F_MOD) $(OBJ_F_FUN) $(OBJ_C_FUN)

#---------------------------------------------------------
#   Default rule to build Fortran modules and functions
#---------------------------------------------------------
#   Note: This doesn't need editing.
#---------------------------------------------------------

# Fortran sources
$(DIR_OBJECT)/%.o: %.f90
	@echo $(FC) $<
	@$(FC) $(OPT_F_COMP) $(PASS_ON) -c -o $@ $<

# C functions
$(DIR_OBJECT)/%.o: %.c
	@echo $(CC) $<
	@$(CC) $(OPT_C_COMP) $(PASS_ON) -c -o $@ $<

# Special rule to build User_Mod.f90.
# Note: If you don't write it as User%Mod.o, makefile will think that it
# is the final target and would not proceed further to build the Program
$(DIR_OBJECT)/User%Mod.o: User%Mod.f90
	@echo $(FC) $<
	@$(FC) $(OPT_F_COMP_USER) $(PASS_ON) -c -o $@ $<

#-----------------------------------
#   Rule to build main program
#-----------------------------------
#   Note: Should not be modified.
#-----------------------------------
$(PROGRAM_FILE): $(OBJ)
	@echo -e Linking "\033[0;32m $(PETSC_LIB) $(METIS_LIB) $(PROGRAM_FILE) \033[0m"
	@$(FC) $(OPT_LINK) -o $(PROGRAM_FILE) $(OBJ) $(METIS_LIB) $(PETSC_LIB)

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
	@rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod \
	       $(DIR_MODULE)/*.f90 $(PROGRAM_FILE)
	@echo done
	@echo -n Removing links to sources in case directory ...
	@rm -f $(DIR_USER)/*.f90
	@echo done
	@echo -n Retreiving default user functions from git repository ...
	@git checkout $(DIR_USER)/*.f90
	@echo done

#-----------------------------------
#   Load options to compile PETSc
#-----------------------------------
FC_OLD := $(FC)
ifneq ($(PETSC_DIR),)
  include ${PETSC_DIR}/lib/petsc/conf/variables
endif

ifneq ($(FC),$(FC_OLD))
$(info #   That might indeed be happening.  Based on your options, or lack     )
$(info #   of any, makefile would expect $(FC_OLD), but PETSc enforces $(FC).  )
$(info #   They might, or might not be compatible, please check.               )
$(info #-----------------------------------------------------------------------)
endif
