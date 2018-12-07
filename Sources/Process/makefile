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
# MPI implementation ("openmpi" or "mpich")
# Warning: change line "FC = mpif90" to your compiler
MPI_TYPE = openmpi

# Link with CGNS ("yes" or "no")
CGNS_HDF5 = no

# Directories for objects and modules. (No need to change.)
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
DIR_MODELS = Turbulence
DIR_USER   = User_Mod
DIR_CASE   = no
VPATH = $(DIR_SHARED):$(DIR_MODELS):.

# Program name (This should hardly change)
PROGRAM_NAME = Process
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with compiler $(FORTRAN)                    )
$(info #-----------------------------------------------------------------------)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=gnu/intel/portland> <DEBUG=yes/no> <CGNS_HDF5=yes/no> )
$(info #        <MPI=yes/no> <MPI_TYPE=openmpi/mpich>                          )
$(info #        <DIR_CASE=relative_path_to_case>                               )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make                       - sequential version with gnu compiler   )
$(info #   make MPI=yes               - parallel version with gnu compiler     )
$(info #   make FORTRAN=intel MPI=yes - parallel version with intel compiler   )
$(info #   make MPI=yes CGNS_HDF5=yes - parallel version with CGNS/HDF5 support)
$(info #   make MPI=yes MPI_TYPE=mpich- parallel version with MPICH/OpenMPI mpi)
$(info #-----------------------------------------------------------------------)

#-------------------------------------------------------------------------------
#   Compiler and linker options
#-------------------------------------------------------------------------------
#   Note: Changes only when support to a new Fortran compiler is added.
#-------------------------------------------------------------------------------

# Fortran == gnu
ifeq ($(FORTRAN), gnu)
  FC = gfortran
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
  FC = ifort
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
  FC = pgfortran
  ifeq ($(DEBUG),yes)
    OPT_COMP = -module $(DIR_MODULE) -r8 -i8 -O0 -g
  else
    OPT_COMP = -module $(DIR_MODULE) -r8 -i8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif

# Use here your command for mpi fortran compiler:
ifeq ($(MPI), yes)
  FC = mpif90
endif

# If compiling with CGNS, add path to include directory
ifeq ($(CGNS_HDF5), yes)
  ifeq ($(MPI), yes)
    ifeq ($(MPI_TYPE), openmpi)
      VALUE = $(shell echo using Cgns library built on OpenMPI MPI implementation)
      DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Openmpi/Include
      CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Openmpi/Lib/libcgns.a \
                    ../Libraries/Hdf5_Par_Openmpi/Lib/libhdf5.a -ldl
    endif
    ifeq ($(MPI_TYPE), mpich)
      VALUE = $(shell echo using Cgns library built on MPICH MPI implementation)
      DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Mpich/Include
      CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Mpich/Lib/libcgns.a \
                    ../Libraries/Hdf5_Par_Mpich/Lib/libhdf5.a -ldl
    endif
  else
    VALUE = $(shell echo using Cgns library built sequentially)
    DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Include
    CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Lib/libcgns.a \
                  ../Libraries/Hdf5_Seq/Lib/libhdf5.a -ldl
  endif
  $(info $(VALUE))
endif
OPT_COMP += $(addprefix -I, $(DIR_CGNS_H))

# METIS and CHACO library
METIS_LIB = ../Libraries/Metis_5.1.0_Linux_64/libmetis.a
CHACO_LIB = ../Libraries/Chaco_2.2_Linux_64/libchaco.a

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
      $(info  $(SRC_USR))

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
SRC_MOD = Const_Mod.f90		\
          Sort_Mod.f90		\
          Tokenizer_Mod.f90	\
          Div_Mod.f90		\
          Grid_Mod.f90		\
          Grid_Level_Mod.f90	\
          Work_Mod.f90		\
          Name_Mod.f90		\
          Material_Mod.f90	\
          Metis_Options_Mod.f90

# Conditional support for MPI communication
ifeq ($(MPI), yes)
  SRC_MOD += Comm_Mod_Par.f90
  COMM_OBJ = Comm_Mod_Par.o
else
  SRC_MOD += Comm_Mod_Seq.f90
  COMM_OBJ = Comm_Mod_Seq.o
endif

# Sources for modules in the local directory
SRC_MOD += Info_Mod.f90		\
           Bulk_Mod.f90		\
           Var_Mod.f90		\
           Matrix_Mod.f90	\
           Vector_Mod.f90	\
           Numerics_Mod.f90	\
           Solver_Mod.f90	\
           Field_Mod.f90	\
           Grad_Mod.f90		\
           Monitor_Mod.f90

# Sources for modules in sub-directories
SRC_MOD += Turbulence_Mod.f90	\
           Les_Mod.f90		\
           Rans_Mod.f90

# Control_Mod should come because it reads Turbulence, Field and Numerics
SRC_MOD += Control_Mod.f90 	\
           Backup_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all functions are obtained by a shell command
SRC_FUN =  $(shell ls -1 *.f90 | grep -v -i _Mod | grep -v -i Save_Cgns | grep -v -i Save_Vtu_Results)
SRC_FUN += $(shell ls -1 $(DIR_MODELS)/*.f90 | xargs -n1 basename | grep -v -i _Mod)

SRC_FUN += Distance.f90		\
           Distance_Squared.f90	\
           Name_File.f90	\
           Approx_Real.f90	\
           Approx_String.f90	\
           Swap_Int.f90		\
           Swap_Real.f90	\
           Swap_Short.f90	\
           To_Upper_Case.f90

#-------------------------------------------
#   Conditional inclusion of CGNS support
#-------------------------------------------
ifeq ($(CGNS_HDF5), no)
  SRC_MOD += Cgns_Mod_Void.f90
  SRC_FUN += Save_Cgns_Results_Void.f90	\
             Save_Cgns_Cells_Void.f90	\
             Save_Vtu_Results.f90
else
  ifeq ($(MPI), yes)
    SRC_MOD += Cgns_Mod_Par.f90
  else
    SRC_MOD += Cgns_Mod_Seq.f90
  endif
  SRC_FUN += Save_Cgns_Results.f90	\
             Save_Cgns_Cells.f90
endif

# Conditional CGNS support is needed for user functions too
ifeq ($(CGNS_HDF5), no)
  SRC_MOD += User_Mod_Vtu.f90
else
  SRC_MOD += User_Mod_Cgns.f90
endif

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
                              $(METIS_LIB) $(CHACO_LIB)

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
