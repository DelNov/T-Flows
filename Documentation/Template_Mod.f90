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
!---------------------------------[Interfaces]---------------------------------!
  interface
    include 'Some_Cool_Header.h90'
    ...
    ..
  end interface
!==============================================================================!

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

  ! Inclusion of member procedures
  contains

  include ’Template_Mod/Allocate.f90’
  ...
  ..

  end module
