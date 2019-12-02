#-------------------------------------------------------------------------------
#
#   Process Makefile
#
#-------------------------------------------------------------------------------

#--------------------------
#   Variable definitions
#--------------------------

# Fortran compiler ("gnu", "intel" or "portland")
FORTRAN = gnu

# Debugging ("yes" or "no")
DEBUG = no

# Compile for MPI ("yes" or "no")
MPI = no

# Link with CGNS ("yes" or "no")
CGNS = no
# MPI implementation used in parallel CGNS lib ("openmpi" or "mpich")
CGNS_MPI = openmpi

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

ifeq ($(MAKECMDGOALS), clean)
$(info #=======================================================================)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=gnu/intel/portland> <DEBUG=yes/no> <MPI=yes/no>       )
$(info #        <FCOMP=gfortran/ifort/pgf90/mpif90/...>                        )
$(info #        <CGNS=yes/no> <CGNS_MPI=openmpi/mpich>                         )
$(info #        <DIR_CASE=full_or_relative_path_to_case>                       )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make                       - sequential version with gnu compiler   )
$(info #   make MPI=yes               - parallel version with gnu compiler     )
$(info #   make MPI=yes FCOMP=compiler- parallel version with custom compiler  )
$(info #   make FORTRAN=intel MPI=yes - parallel version with intel compiler   )
$(info #   make MPI=yes CGNS=yes      - parallel version with CGNS/HDF5 support)
$(info #-----------------------------------------------------------------------)
endif
#-------------------------------------------------------------------------------
#   Compiler and linker options
#-------------------------------------------------------------------------------
#   Note: Changes only when support to a new Fortran compiler is added.
#-------------------------------------------------------------------------------

# Fortran == gnu
ifeq ($(FORTRAN), gnu)
  ifeq ($(FCOMP),)
    FC = gfortran
  else
    FC = $(FCOMP)
  endif
  ifeq ($(DEBUG),yes)
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -fdefault-integer-8 -O0 -g \
               -ffixed-line-length-80 -fcheck=all -fwhole-file -fbacktrace  \
               -ffpe-trap=invalid,zero,overflow                             \
               -fimplicit-none -finit-real=nan                              \
               -Wall -Wextra -Waliasing -Wampersand -Warray-temporaries     \
               -Wc-binding-type -Wcharacter-truncation -Wline-truncation    \
               -Wconversion -Wconversion-extra -Wintrinsics-std             \
               -Wreal-q-constant -Wsurprising -Wtabs -Wunderflow            \
               -Wintrinsic-shadow -Wunused-parameter -Walign-commons        \
               -Wfunction-elimination -Wrealloc-lhs -Wrealloc-lhs-all       \
               -Wcompare-reals -Wtarget-lifetime -Wpedantic -fstack-check
  else
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -fdefault-integer-8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif

# Fortran == intel
ifeq ($(FORTRAN), intel)
  ifeq ($(FCOMP),)
    FC = ifort
  else
    FC = $(FCOMP)
  endif
  ifeq ($(DEBUG),yes)
    OPT_COMP = -module $(DIR_MODULE) -r8 -i8 -O0 -g -warn all -check all \
               -debug all -fpe-all=0 -traceback
  else
    OPT_COMP = -module $(DIR_MODULE) -r8 -i8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif

# Fortran == portland
ifeq ($(FORTRAN), portland)
  ifeq ($(FCOMP),)
    FC = pgfortran
  else
    FC = $(FCOMP)
  endif
  ifeq ($(DEBUG),yes)
    OPT_COMP = -module $(DIR_MODULE) -r8 -i8 -O0 -g
  else
    OPT_COMP = -module $(DIR_MODULE) -r8 -i8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif

# Use here your command for mpi fortran compiler:
ifeq ($(MPI), yes)
  ifeq ($(FCOMP),)
    FC = mpif90
  else
    FC = $(FCOMP)
  endif
endif

# If not "make clean"
ifneq ($(MAKECMDGOALS), clean)
$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with $(FC) compiler                         )
$(info #-----------------------------------------------------------------------)
endif

# If compiling with CGNS, add path to include directory
ifeq ($(CGNS), yes)
  ifeq ($(MPI), yes)
    ifeq ($(CGNS_MPI), openmpi)
      VALUE = $(shell echo Using CGNS library built on OpenMPI MPI implementation)
      DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Openmpi/Include
      CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Openmpi/Lib/libcgns.a \
                    ../Libraries/Hdf5_Par_Openmpi/Lib/libhdf5.a -ldl
    endif
    ifeq ($(CGNS_MPI), mpich)
      VALUE = $(shell echo Using CGNS library built on MPICH MPI implementation)
      DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Mpich/Include
      CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Mpich/Lib/libcgns.a \
                    ../Libraries/Hdf5_Par_Mpich/Lib/libhdf5.a -ldl
    endif
  else
    VALUE = $(shell echo Using CGNS library built sequentially)
    DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Include
    CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Lib/libcgns.a \
                  ../Libraries/Hdf5_Seq/Lib/libhdf5.a -ldl
  endif
  $(info $(VALUE))
endif
OPT_COMP += $(addprefix -I, $(DIR_CGNS_H))

# METIS library
METIS_LIB = ../Libraries/Metis_5.1.0_Linux_64/libmetis.a

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

# Conditional support for MPI communication
ifeq ($(MPI), yes)
  SRC_MOD  = Comm_Mod_Par.f90
else
  SRC_MOD  = Comm_Mod_Seq.f90
endif

# Modules in shared directory
SRC_MOD += Const_Mod.f90	\
           Cpu_Timer_Mod.f90	\
           Sort_Mod.f90		\
           Math_Mod.f90		\
           File_Mod.f90		\
           Div_Mod.f90		\
           Grid_Mod.f90		\
           Grid_Level_Mod.f90	\
           Work_Mod.f90		\
           Material_Mod.f90	\
           Metis_Options_Mod.f90

# Sources for modules in the local directory
SRC_MOD += Control_Mod.f90	\
           Info_Mod.f90		\
           Bulk_Mod.f90		\
           Var_Mod.f90		\
           Face_Mod.f90		\
           Matrix_Mod.f90	\
           Vector_Mod.f90	\
           Solver_Mod.f90	\
           Numerics_Mod.f90	\
           Field_Mod.f90	\
           Eddies_Mod.f90	\
           Monitor_Mod.f90	\
           Turb_Mod.f90		\
           Swarm_Mod.f90	\
           Surf_Mod.f90		\
           Multiphase_Mod.f90

# User_Mod should come before Save_Results_Mod
SRC_MOD += User_Mod.f90

# Save_Results_Mod - conditional VTU/CGNS support
ifeq ($(CGNS), no)
  SRC_MOD += Save_Grid_Mod_Vtu.f90	\
             Save_Results_Mod_Vtu.f90
else
  ifeq ($(MPI), yes)
    SRC_MOD += Cgns_Mod_Par.f90
  else
    SRC_MOD += Cgns_Mod_Seq.f90
  endif
  SRC_MOD += Save_Grid_Mod_Cgns.f90	\
             Save_Results_Mod_Cgns.f90
endif

# Backup_Mod should come last because it uses all other modules
SRC_MOD += Backup_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all functions in the local ...
# ... directory are obtained by a shell command
SRC_FUN =  $(shell ls -1 *.f90 | grep -v -i _Mod)

# Non-member functions still residing in Turb_Mod
SRC_FUN += Tau_Wall_Low_Re.f90		\
           Tau_Wall_Rough_Walls.f90	\
           Y_Plus_Low_Re.f90		\
           Y_Plus_Rough_Walls.f90	\
           U_Plus_Log_Law.f90		\
           U_Plus_Rough_Walls.f90	\
           Time_And_Length_Scale.f90	\
           Find_Nearest_Wall_Cell.f90	\
           Roughness_Coefficient.f90

SRC_FUN += Swap_Int.f90		\
           Swap_Real.f90	\
           Swap_Short.f90	\
           To_Upper_Case.f90

#----------------------------------------------------------------------
#   List of objects generated from the list of modules and functions
#----------------------------------------------------------------------
#   Note: This doesn't need editing.
#----------------------------------------------------------------------
OBJ_MOD = $(SRC_MOD:%.f90=$(DIR_OBJECT)/%.o)
OBJ_FUN = $(SRC_FUN:%.f90=$(DIR_OBJECT)/%.o)
OBJ = $(OBJ_MOD) $(OBJ_FUN)

#-------------------------------------------------------
#   List of modules currently used for target "clean"
#-------------------------------------------------------
#   Note: This doesn't need editing.
#-------------------------------------------------------
SRC_MOD_LOW = $(shell echo $(SRC_MOD) | tr A-Z a-z)
MOD = $(SRC_MOD_LOW:%.f90=$(DIR_MODULE)/%.mod)

#---------------------------------------------------------
#   Default rule to build Fortran modules and functions
#---------------------------------------------------------
#   Note: This doesn't need editing.
#---------------------------------------------------------

# Modules
$(DIR_OBJECT)/%.o: %.f90 %/*.f90
	@echo FC $<
	@$(FC) $(OPT_COMP) -c -o $@ $<

# Functions
$(DIR_OBJECT)/%.o: %.f90
	@echo FC $<
	@$(FC) $(OPT_COMP) -c -o $@ $<

#-----------------------------------
#   Rule to build main program
#-----------------------------------
#   Note: Should not be modified.
#-----------------------------------
$(PROGRAM_FILE): $(OBJ)
	@echo Linking "\033[0;32m $(PROGRAM_FILE) \033[0m"
	@$(FC) $(OPT_LINK) -o $(PROGRAM_FILE) $(OBJ) $(CGNS_LIB)  \
                              $(METIS_LIB)

#--------------------------------------------------------------------
#   Explicit dependencies for modules
#--------------------------------------------------------------------
#   These are automatically generated by:
#   Sources/Utilities/create_external_dependencies_for_makefile.sh
#--------------------------------------------------------------------
include makefile_explicit_dependencies

#---------------------
#   Explicit target.
#---------------------
clean:
	@echo -n Removing objects and modules ... 
	@rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
	@echo done
	@echo -n Removing links to sources in case directory ... 
	@rm -f $(DIR_USER)/*.f90
	@echo done
	@echo -n Retreiving default user functions from git repository ... 
	@git checkout $(DIR_USER)/*.f90
	@echo done
