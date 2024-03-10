#include "Unused.h90"

!==============================================================================!
  module Gpu_Mod
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Native_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------!
  !              !
  !   GPU type   !
  !              !
  !--------------!
  type Gpu_Type

    real :: gb_used = 0.0

    contains

      ! Procedures to copy vectors (and matrices) to device
      procedure :: Vector_Int_Copy_To_Device
      procedure :: Vector_Int_Create_On_Device
      procedure :: Vector_Int_Destroy_On_Device
      procedure :: Vector_Real_Copy_To_Device
      procedure :: Vector_Real_Create_On_Device
      procedure :: Vector_Real_Destroy_On_Device
      procedure :: Vector_Update_Device
      procedure :: Vector_Update_Host
      procedure :: Matrix_Int_Copy_To_Device
      procedure :: Matrix_Int_Destroy_On_Device
      procedure :: Matrix_Real_Copy_To_Device
      procedure :: Matrix_Real_Destroy_On_Device

      ! Procedures to copy sparse matrices to device
      procedure :: Sparse_Copy_To_Device
      procedure :: Sparse_Destroy_On_Device

      ! Procedures to copy native solver to device
      procedure :: Native_Transfer_To_Device
      procedure :: Native_Destroy_On_Device

  end type

  type(Gpu_Type) :: Gpu

  contains

    ! Procedures to copy vectors and matrices to device
#   include "Gpu_Mod/Vector/Int_Copy_To_Device.f90"
#   include "Gpu_Mod/Vector/Int_Create_On_Device.f90"
#   include "Gpu_Mod/Vector/Int_Destroy_On_Device.f90"
#   include "Gpu_Mod/Vector/Real_Copy_To_Device.f90"
#   include "Gpu_Mod/Vector/Real_Create_On_Device.f90"
#   include "Gpu_Mod/Vector/Real_Destroy_On_Device.f90"
#   include "Gpu_Mod/Vector/Update_Device.f90"
#   include "Gpu_Mod/Vector/Update_Host.f90"
#   include "Gpu_Mod/Matrix/Int_Copy_To_Device.f90"
#   include "Gpu_Mod/Matrix/Int_Destroy_On_Device.f90"
#   include "Gpu_Mod/Matrix/Real_Copy_To_Device.f90"
#   include "Gpu_Mod/Matrix/Real_Destroy_On_Device.f90"

    ! Procedures to copy vectors to device
#   include "Gpu_Mod/Sparse/Copy_To_Device.f90"
#   include "Gpu_Mod/Sparse/Destroy_On_Device.f90"

    ! Procedures to copy native solver to device
#   include "Gpu_Mod/Native/Transfer_To_Device.f90"
#   include "Gpu_Mod/Native/Destroy_On_Device.f90"

  end module
