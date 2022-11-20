#===============================================================================
#
#   Divide Makefile
#
#-------------------------------------------------------------------------------

#-------------------------------------
#   Default options for compilation
#-------------------------------------
FORTRAN ?= gnu
DEBUG   ?= no
ASSERT  ?= yes
REAL    ?= double
MPI     ?= no
OPENMP  ?= no

# Directories for objects and modules. (No need to change.)
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
VPATH = $(DIR_SHARED):.

# Program name (This should hardly change)
PROGRAM_NAME = Divide
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

#----------------------------------------------------------
#   Initialize variables to pass on to preprocessor with
#   variable T_FLOWS_COMPILATION which, if set to one,
#   means unit is compiled from T-Flow
#   It later adds T_FLOWS_MPI and T_FLOWS_PETSC
#----------------------------------------------------------
PASS_ON = -DT_FLOWS_COMPILATION=1
ifeq ($(ASSERT), yes)
  PASS_ON += -DT_FLOWS_ASSERT=1
else
  PASS_ON += -DT_FLOWS_ASSERT=0
endif

$(info #=======================================================================)
$(info # Compiling $(PROGRAM_NAME) with compiler $(FORTRAN)                    )
$(info #-----------------------------------------------------------------------)
$(info # Usage:                                                                )
$(info #   make <FORTRAN=gnu/intel/nvidia> <DEBUG=no/yes> <ASSERT=yes/no>      )
$(info #        <FCOMP=gfortran/ifort/nvfortran/mpif90/mpiifort/...>           )
$(info #        <REAL=double/single>                                           )
$(info #                                                                       )
$(info # Note: The first item, for each of the options above, is the default.  )
$(info #                                                                       )
$(info # Examples:                                                             )
$(info #   make              - compile with gnu compiler                       )
$(info #   make FORTAN=intel - compile with intel compiler                     )
$(info #   make DEBUG=yes    - compile with gnu compiler in debug mode         )
$(info #   make ASSERT=no    - compile without assert with gnu compiler        )
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
                 -check all -debug all -fpe-all=0 -traceback -cpp
  else
    OPT_F_COMP = -module $(DIR_MODULE) -O3 -cpp
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -r8
  endif
  ifeq ($(OPENMP), yes)
    OPT_F_COMP += -qopenmp
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
  ifeq ($(OPENMP), yes)
    OPT_F_COMP += -mp
  endif

# Fortran == gnu
else
  ifeq ($(FCOMP),)
    FC = gfortran
  else
    FC = $(FCOMP)
  endif
  ifeq ($(DEBUG), yes)
    OPT_F_COMP = -J $(DIR_MODULE) -O0 -g -cpp                                  \
                 -ffree-line-length-none -fcheck=all -fwhole-file -fbacktrace  \
                 -ffpe-trap=invalid,zero,overflow -fimplicit-none              \
                 -Wall -Wextra -Waliasing -Wampersand -Warray-temporaries      \
                 -Wc-binding-type -Wcharacter-truncation -Wline-truncation     \
                 -Wconversion -finit-real=nan    -Wintrinsics-std              \
                 -Wreal-q-constant -Wsurprising -Wtabs -Wunderflow             \
                 -Wintrinsic-shadow -Wunused-parameter -Walign-commons         \
                 -Wfunction-elimination -Wrealloc-lhs -Wrealloc-lhs-all        \
                 -Wtarget-lifetime -Wpedantic -fstack-check
  else
    OPT_F_COMP = -J $(DIR_MODULE) -O3 -cpp -ffree-line-length-none
  endif
  ifeq ($(REAL), double)
    OPT_F_COMP += -fdefault-real-8
  endif
endif

#-----------------------------------------------------------
#   Set linker options to be the same as compiler options
#-----------------------------------------------------------
OPT_LINK = $(OPT_F_COMP)

#-----------------------
#   MPI related stuff
#-----------------------

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
#   Note: Modules written in lower case
#         letters are candidates for deletion.
#------------------------------------------------------

#-------------
#   Modules
#-------------

# Modules in shared directories
SRC_F_MOD = Const_Mod.f90      \
            String_Mod.f90     \
            Tokenizer_Mod.f90  \
            Comm_Mod.f90       \
            Assert_Mod.f90     \
            Message_Mod.f90    \
            Profiler_Mod.f90   \
            Math_Mod.f90       \
            File_Mod.f90       \
            Gen_Mod.f90        \
            Vtk_Mod.f90        \
            Metis_Mod.f90      \
            Boundary_Mod.f90   \
            Swap_Mod.f90       \
            Sort_Mod.f90       \
            Grid_Mod.f90

# Modules in the local directory
SRC_F_MOD += Divide_Mod.f90

#---------------
#   Functions
#---------------

# Sources for all shared functions
SRC_F_FUN = Adjust_Dim.f90        \
            Adjust_First_Dim.f90

# Local functions are added one by one ...
SRC_F_FUN += Main_Div.f90

#----------------------------------------------------------------------
#   List of objects generated from the list of modules and functions
#----------------------------------------------------------------------
#   Note: This doesn't need editing.
#----------------------------------------------------------------------
OBJ_F_MOD = $(SRC_F_MOD:%.f90=$(DIR_OBJECT)/%.o)
OBJ_F_FUN = $(SRC_F_FUN:%.f90=$(DIR_OBJECT)/%.o)
OBJ = $(OBJ_F_MOD) $(OBJ_F_FUN)

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
	@$(FC) $(OPT_F_COMP) $(PASS_ON) -c -o $@ $<

# Fortran functions
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
