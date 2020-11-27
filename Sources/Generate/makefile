#===============================================================================
#
#   Generate Makefile
#
#-------------------------------------------------------------------------------

#--------------------------
#   Variable definitions
#--------------------------

# Fortran compiler ("gnu", "intel" or "portland")
FORTRAN = gnu

# Debugging ("yes" or "no")
DEBUG = no

# Link with CGNS ("yes" or "no")
CGNS = no

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
$(info #   make <FORTRAN=gnu/intel/portland> <DEBUG=yes/no> <CGNS=yes/no>      )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make              - compile with gnu compiler                       )
$(info #   make FORTAN=intel - compile with intel compiler                     )
$(info #   make DEBUG=yes    - compile with gnu compiler in debug mode         )
$(info #   make CGNS=yes     - compile with gnu compiler and CGNS support      )
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

# If compiling with CGNS, add path to include directory
ifeq ($(CGNS), yes)
  DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Include
  CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Lib/libcgns.a \
                ../Libraries/Hdf5_Seq/Lib/libhdf5.a -ldl
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

#-------------
#   Modules
#-------------

# Modules in shared directories
SRC_MOD = Const_Mod.f90			\
          Comm_Mod_Seq.f90		\
          Math_Mod.f90			\
          File_Mod.f90			\
          Gen_Mod.f90			\
          Metis_Options_Mod.f90		\
          Bnd_Cond_Mod.f90		\
          Sort_Mod.f90			\
          Grid_Mod.f90

# Modules in the local directory
SRC_MOD += Point_Mod.f90		\
           Line_Mod.f90			\
           Block_Mod.f90		\
           Region_Mod.f90		\
           Smooths_Mod.f90		\
           Refines_Mod.f90		\
           Domain_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all shared functions
SRC_FUN = Probe_1d_Nodes.f90	\
          Probe_2d.f90		\
          Swap_Int.f90		\
          Swap_Real.f90		\
          Swap_Short.f90	\
          To_Upper_Case.f90

SRC_FUN += Calculate_Geometry.f90	\
           Load_Dom.f90			\
           Logo_Gen.f90			\
           Main_Gen.f90			\
           Print_Grid_Statistics.f90

#-------------------------------------------
#   Conditional inclusion of CGNS support
#-------------------------------------------
ifeq ($(CGNS), no)
  SRC_MOD += Save_Grid_Mod_Vtu.f90
else
  SRC_MOD += Cgns_Mod_Seq.f90
  SRC_MOD += Save_Grid_Mod_Cgns.f90
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
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
