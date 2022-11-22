!==============================================================================!
  module Stl_Mod
!------------------------------------------------------------------------------!
!   Module to describe geometry in STL format                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
!==============================================================================!

  !--------------!
  !   Stl type   !
  !--------------!
  type Stl_Type
    character(SL)     :: name                    ! file name which defines it
    integer           :: n_facets                ! number of facets
    real, allocatable :: x(:,:), y(:,:), z(:,:)  ! three vertex's coordinates
    real, allocatable :: nx(:),  ny(:),  nz(:)   ! normal vectors

    contains
      procedure :: Create_From_File

  end type

  contains
#   include "Stl_Mod/Create_From_File.f90"

  end module
