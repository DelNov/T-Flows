#-------------------------------------------------------------------------------
#
#   Process Makefile
#
#-------------------------------------------------------------------------------

#-------------------------------------
#   Default options for compilation
#-------------------------------------
FORTRAN ?= gnu
DEBUG   ?= no
REAL    ?= double
MPI     ?= no

# Directories for objects and modules. (No need to change.)
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
DIR_TURBUL = Turb_Mod
DIR_SWARM  = Swarm_Mod
DIR_USER   = User_Mod
DIR_CASE   = no
VPATH = $(DIR_SHARED):$(DIR_TURBUL):$(DIR_SWARM):.

# Program name (This should hardly change)
PROGRAM_NAME = Process
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)


$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with $(FORTRAN) compiler                    )
$(info #-----------------------------------------------------------------------)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=gnu/intel/nvidia> <DEBUG=no/yes> <MPI=no/yes>         )
$(info #        <FCOMP=gfortran/ifort/nvfortran/mpif90/...>                    )
$(info #        <DIR_CASE=full_or_relative_path_to_case>                       )
$(info #        <REAL=double/single>                                           )
$(info #                                                                       )
$(info # Notes:                                                                )
$(info #   1. The first item, for each of the options above, is the default.   )
$(info #   2. If PETSC_DIR environment variable is defined, the code willl     )
$(info #      make an attempt to compile with PETSc libraries.                 )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make                        - sequential version with gnu compiler  )
$(info #   make MPI=yes                - MPI version with gnu compiler         )
$(info #   make MPI=yes FCOMP=compiler - MPI version with custom compiler      )
$(info #   make FORTRAN=intel MPI=yes  - MPI version with intel compiler       )
$(info #-----------------------------------------------------------------------)

#-------------------
#   METIS library
#-------------------
METIS_LIB = ../Libraries/Metis_5.1.0_Linux_64/libmetis_i32_r64.a
ifeq ($(REAL), single)
  METIS_LIB = ../Libraries/Metis_5.1.0_Linux_64/libmetis_i32_r32.a
endif

#-------------------------------------------------------------------------
#   Compiler and linker options
#-------------------------------------------------------------------------
#   Note: Changes only when support to a new Fortran compiler is added.
#-------------------------------------------------------------------------

# Fortran == intel
ifeq ($(FORTRAN), intel)
  ifeq ($(FCOMP),)
    FC = ifort
  else
    FC = $(FCOMP)
  endif
  ifeq ($(DEBUG), yes)
    OPT_F_COMP = -module $(DIR_MODULE) -O0 -g -warn -warn all,noexternal  \
                 -check all -debug all -fpe-all=0 -traceback
  else
    OPT_F_COMP = -module $(DIR_MODULE) -O3
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -r8
  endif

# Fortran == nvidia
else ifeq ($(FORTRAN), nvidia)
  ifeq ($(FCOMP),)
    FC = nvfortran
  else
    FC = $(FCOMP)
  endif
  ifeq ($(DEBUG), yes)
    OPT_F_COMP = -module $(DIR_MODULE) -O0 -g
  else
    OPT_F_COMP = -module $(DIR_MODULE) -O3
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -r8
  endif

# Fortran == gnu
else
  ifeq ($(FCOMP),)
    FC = gfortran
  else
    FC = $(FCOMP)
  endif
  ifeq ($(DEBUG), yes)
    OPT_F_COMP = -J $(DIR_MODULE) -O0 -g -cpp                                 \
                 -ffixed-line-length-80 -fcheck=all -fwhole-file -fbacktrace  \
                 -ffpe-trap=invalid,zero,overflow -fimplicit-none             \
                 -Wall -Wextra -Waliasing -Wampersand -Warray-temporaries     \
                 -Wc-binding-type -Wcharacter-truncation -Wline-truncation    \
                 -Wconversion -finit-real=nan    -Wintrinsics-std             \
                 -Wreal-q-constant -Wsurprising -Wtabs -Wunderflow            \
                 -Wintrinsic-shadow -Wunused-parameter -Walign-commons        \
                 -Wfunction-elimination -Wrealloc-lhs -Wrealloc-lhs-all       \
                 -Wtarget-lifetime -Wpedantic -fstack-check
  else
    OPT_F_COMP = -J $(DIR_MODULE) -O3 -cpp
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -fdefault-real-8
  endif
endif

# Set initial options for C compiler - this is probably not needed
OPT_C_COMP = -cpp

#----------------------------------------------
#   Additional options to compile with PETSc
#----------------------------------------------
ifneq ($(PETSC_DIR),)
  OPT_F_COMP += -I ${PETSC_DIR}/include
  OPT_F_COMP += -I ${PETSC_DIR}/${PETSC_ARCH}/include
  OPT_C_COMP += -I ${PETSC_DIR}/include
  OPT_C_COMP += -I ${PETSC_DIR}/${PETSC_ARCH}/include
endif

#-----------------------------------------------------------
#   Set linker options to be the same as compiler options
#-----------------------------------------------------------
OPT_LINK = $(OPT_F_COMP)

# Use here your command for mpi fortran compiler:
ifeq ($(MPI), yes)
  ifeq ($(FCOMP),)
    FC = mpif90
  else
    FC = $(FCOMP)
  endif
endif

#------------------------------------------------------
#   List of sources for modules and functions
#------------------------------------------------------
#   Modules' order must obey their dependency
#   This list should therefore be written "by hand".
#   Note: Modules written in lower case
#         letters are candidates for deletion.
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

# Const is the lowest level module
SRC_F_MOD = Const_Mod.f90

# Conditional support for MPI communication
ifeq ($(MPI), yes)
  SRC_F_MOD += Comm_Mod_Par.f90
else
  SRC_F_MOD += Comm_Mod_Seq.f90
endif

# Modules in shared directory
SRC_F_MOD += Cpu_Timer_Mod.f90      \
             Math_Mod.f90           \
             Swap_Mod.f90           \
             Sort_Mod.f90           \
             File_Mod.f90           \
             Control_Mod.f90        \
             Vtk_Mod.f90            \
             Boundary_Mod.f90       \
             Metis_Options_Mod.f90  \
             Grid_Mod.f90           \
             Work_Mod.f90

# Sources for modules in the local directory
SRC_F_MOD += Point_Mod.f90       \
             Bulk_Mod.f90        \
             Var_Mod.f90         \
             Face_Mod.f90        \
             Vert_Mod.f90        \
             Elem_Mod.f90        \
             Side_Mod.f90        \
             Matrix_Mod.f90      \
             Vector_Mod.f90      \
             Native_Mod.f90
ifneq ($(PETSC_DIR),)
  SRC_F_MOD += Petsc_Mod_True.f90
  SRC_C_FUN = C_Petsc_Interface_True.c
else
  SRC_F_MOD += Petsc_Mod_Fake.f90
  SRC_C_FUN = C_Petsc_Interface_Fake.c
endif
SRC_F_MOD += Solver_Mod.f90      \
             Numerics_Mod.f90    \
             Field_Mod.f90       \
             Front_Mod.f90       \
             Surf_Mod.f90        \
             Info_Mod.f90        \
             Eddies_Mod.f90      \
             Monitor_Mod.f90     \
             Turb_Mod.f90        \
             Vof_Mod.f90         \
             Particle_Mod.f90    \
             Swarm_Mod.f90       \
             Interface_Mod.f90

# User_Mod should come before Results_Mod
SRC_F_MOD += User_Mod.f90

# Backup_Mod should come close to the end because it uses all other modules
# Note the conditional support for MPI communication
ifeq ($(MPI), yes)
  SRC_F_MOD += Backup_Mod_Par.f90
else
  SRC_F_MOD += Backup_Mod_Seq.f90
endif

# Results_Mod must be last.
SRC_F_MOD += Results_Mod.f90       \
             Read_Controls_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all functions in the local directory
SRC_F_FUN = Balance_Volume.f90                 \
            Calculate_Shear_And_Vorticity.f90  \
            Calculate_Vorticity.f90            \
            Compute_Energy.f90                 \
            Compute_Momentum_Explicit.f90      \
            Compute_Momentum.f90               \
            Compute_Pressure.f90               \
            Compute_Scalar.f90                 \
            Convective_Outflow.f90             \
            Correct_Velocity.f90               \
            Initialize_Variables.f90           \
            Key_Ind.f90                        \
            Logo_Pro.f90                       \
            Main_Pro.f90                       \
            Piso_Algorithm.f90                 \
            Rhie_And_Chow.f90                  \
            Update_Boundary_Values.f90

SRC_F_FUN += 

#----------------------------------------------------------------------
#   List of objects generated from the list of modules and functions
#----------------------------------------------------------------------
#   Note: This doesn't need editing.
#----------------------------------------------------------------------
OBJ_F_MOD = $(SRC_F_MOD:%.f90=$(DIR_OBJECT)/%.o)
OBJ_F_FUN = $(SRC_F_FUN:%.f90=$(DIR_OBJECT)/%.o)
OBJ_C_FUN = $(SRC_C_FUN:%.c=$(DIR_OBJECT)/%.o)
OBJ = $(OBJ_F_MOD) $(OBJ_F_FUN) $(OBJ_C_FUN)

#-------------------------------------------------------
#   List of modules currently used for target "clean"
#-------------------------------------------------------
#   Note: This doesn't need editing.
#-------------------------------------------------------
SRC_F_MOD_LOW = $(shell echo $(SRC_F_MOD) | tr A-Z a-z)
MOD = $(SRC_F_MOD_LOW:%.f90=$(DIR_MODULE)/%.mod)

#---------------------------------------------------------
#   Default rule to build Fortran modules and functions
#---------------------------------------------------------
#   Note: This doesn't need editing.
#---------------------------------------------------------

# Fortran modules
$(DIR_OBJECT)/%.o: %.f90 %/*.f90
	@echo $(FC) $<
	@$(FC) $(OPT_F_COMP) -c -o $@ $<

# Fortran functions
$(DIR_OBJECT)/%.o: %.f90
	@echo $(FC) $<
	@$(FC) $(OPT_F_COMP) -c -o $@ $<

# C functions
$(DIR_OBJECT)/%.o: %.c
	@echo $(CC) $<
	@$(CC) $(OPT_C_COMP) -c -o $@ $<

#-----------------------------------
#   Rule to build main program
#-----------------------------------
#   Note: Should not be modified.
#-----------------------------------
$(PROGRAM_FILE): $(OBJ)
	@echo Linking "\033[0;32m $(PETSC_LIB) $(METIS_LIB) $(PROGRAM_FILE) \033[0m"
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
ifneq ($(PETSC_DIR),)
  include ${PETSC_DIR}/lib/petsc/conf/variables
endif
