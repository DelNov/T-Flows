#include "../../Shared/Assert.h90"  ! path to Shared directory may be different
#include "../../Shared/Browse.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Template_Mod
!------------------------------------------------------------------------------!
!   Template for a module                                                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use This_Mod
  use That_Mod
  ...
  ..
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Parameters used inside the module
  integer, parameter :: MAX_ARRAY_SIZE = 1024

  !-------------------!
  !   Template type   !
  !-------------------!
  type Template_Type

  ! Declaration of member data
  ...
  ..

  ! Declaration of member procedures
  contains
    procedure          :: Allocate
    procedure          :: Deploy
    procedure          :: Destroy
    procedure, private :: Use_Only_From_Here

  end type

  ! Declaration of singleton objects, if any, comes next.  In this
  ! example, object Template (of Template_Type) will be available
  ! to other modules which use Template_Mod.  Typical examples are
  ! File and Line defined in File_Mod
  type(Template_Type) :: Template

  ! Inclusion of member procedures
  contains

#  include "Template_Mod/Allocate.f90"
#  include "Template_Mod/Deploy.f90"
#  include "Template_Mod/Destroy.f90"
  ...
  ..

  end module
