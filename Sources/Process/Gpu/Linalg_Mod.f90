# include "Unused.h90"

!==============================================================================!
  module Linalg_Mod
!----------------------------------[Modules]-----------------------------------!
  use Sparse_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !                  !
  !   Compute type   !
  !                  !
  !------------------!
  type Linalg_Type

    contains
      procedure :: Mat_X_Vec
      procedure :: Sca_O_Dia
      procedure :: Set_Singular
      procedure :: Sys_Normalize
      procedure :: Sys_Restore
      procedure :: Vec_Copy
      procedure :: Vec_D_Vec
      procedure :: Vec_M_Vec_X_Vec
      procedure :: Vec_P_Sca_X_Vec
      procedure :: Vec_P_Vec_X_Vec
      procedure :: Vec_X_Vec

  end type

  !-----------------------------------!
  !   Singletone type global object   !
  !-----------------------------------!
  type(Linalg_Type) :: Linalg

  contains
#   include "Linalg_Mod/Mat_X_Vec.f90"
#   include "Linalg_Mod/Sca_O_Dia.f90"
#   include "Linalg_Mod/Set_Singular.f90"
#   include "Linalg_Mod/Sys_Normalize.f90"
#   include "Linalg_Mod/Sys_Restore.f90"
#   include "Linalg_Mod/Vec_Copy.f90"
#   include "Linalg_Mod/Vec_D_Vec.f90"
#   include "Linalg_Mod/Vec_M_Vec_X_Vec.f90"
#   include "Linalg_Mod/Vec_P_Sca_X_Vec.f90"
#   include "Linalg_Mod/Vec_P_Vec_X_Vec.f90"
#   include "Linalg_Mod/Vec_X_Vec.f90"

  end module
