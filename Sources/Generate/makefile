#===============================================================================
#
#   Generate Makefile
#
#-------------------------------------------------------------------------------

#--------------------------
#   Variable definitions
#--------------------------

# Fortran compiler ("gnu", "intel" or "nvidia")
FORTRAN = gnu

# Debugging ("yes" or "no")
DEBUG = no

# Directories for objects and modules. (No need to change.)
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
VPATH = $(DIR_SHARED):.

# Program name (This should hardly change)
PROGRAM_NAME = Generate
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with compiler $(FORTRAN)                    )
$(info #-----------------------------------------------------------------------)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=gnu/intel/nvidia> <DEBUG=yes/no>                      )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make              - compile with gnu compiler                       )
$(info #   make FORTAN=intel - compile with intel compiler                     )
$(info #   make DEBUG=yes    - compile with gnu compiler in debug mode         )
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
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -O0 -g                     \
               -ffixed-line-length-80 -fcheck=all -fwhole-file -fbacktrace  \
               -ffpe-trap=invalid,zero,overflow -fimplicit-none             \
               -Wall -Wextra -Waliasing -Wampersand -Warray-temporaries     \
               -Wc-binding-type -Wcharacter-truncation -Wline-truncation    \
               -Wconversion -finitreal=nan     -Wintrinsics-std             \
               -Wreal-q-constant -Wsurprising -Wtabs -Wunderflow            \
               -Wintrinsic-shadow -Wunused-parameter -Walign-commons        \
               -Wfunction-elimination -Wrealloc-lhs -Wrealloc-lhs-all       \
               -Wcompare-reals -Wtarget-lifetime -Wpedantic -fstack-check
  else
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif

# Fortran == intel
ifeq ($(FORTRAN), intel)
  FC = ifort
  ifeq ($(DEBUG), yes)
    OPT_COMP = -module $(DIR_MODULE) -r8 -O0 -g -warn all -check all \
               -debug all -fpe-all=0 -traceback
  else
    OPT_COMP = -module $(DIR_MODULE) -r8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif

# Fortran == nvidia
ifeq ($(FORTRAN), nvidia)
  FC = nvfortran
  ifeq ($(DEBUG), yes)
    OPT_COMP = -module $(DIR_MODULE) -r8 -O0 -g
  else
    OPT_COMP = -module $(DIR_MODULE) -r8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif

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

#-------------
#   Modules
#-------------

# Modules in shared directories
SRC_MOD = Const_Mod.f90          \
          Comm_Mod_Seq.f90       \
          Math_Mod.f90           \
          File_Mod.f90           \
          Gen_Mod.f90            \
          Vtk_Mod.f90            \
          Metis_Options_Mod.f90  \
          Boundary_Mod.f90       \
          Swap_Mod.f90           \
          Sort_Mod.f90           \
          Grid_Mod.f90

# Modules in the local directory
SRC_MOD += Point_Mod.f90    \
           Line_Mod.f90     \
           Block_Mod.f90    \
           Region_Mod.f90   \
           Smooths_Mod.f90  \
           Refines_Mod.f90  \
           Domain_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all shared functions
SRC_FUN = Probe_1d_Nodes.f90  \
          Probe_2d.f90        \
          To_Upper_Case.f90

# Local functions are added one by one ...
SRC_FUN += Calculate_Geometry.f90  \
           Load_Dom.f90            \
           Logo_Gen.f90            \
           Main_Gen.f90            \
           Print_Grid_Statistics.f90

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
	@$(FC) $(OPT_LINK) -o $(PROGRAM_FILE) $(OBJ) $(METIS_LIB)

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
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
