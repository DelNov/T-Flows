#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Macros.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Gpu_Mod
!----------------------------------[Modules]-----------------------------------!
!>  Handles GPU data transfer and management for T-Flows.
!>
!>  This module provides the infrastructure for handling all GPU-related
!>  data transfer within T-Flows. It manages the memory on the GPU, based on
!>  allocation, copying, and destruction of various T-Flows' data structures
!>  such as vectors, matrices, fields, and grid data.
!>
!>  Each procedure within the module is dedicated to a specific type of
!>  data or operation, ensuring efficient interaction with the GPU.
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
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
      procedure :: Sparse_Con_Copy_To_Device
      procedure :: Sparse_Con_Destroy_On_Device
      procedure :: Sparse_Val_Copy_To_Device
      procedure :: Sparse_Val_Destroy_On_Device

      ! Procedures to copy native solver to device
      procedure :: Native_Copy_To_Device
      procedure :: Native_Destroy_On_Device

      ! Procedures to copy field to device
      procedure :: Field_Copy_To_Device
      procedure :: Field_Destroy_On_Device
      procedure :: Field_Update_Host

      ! Procedures to copy grid to device
      procedure :: Grid_Copy_To_Device
      procedure :: Grid_Destroy_On_Device
      procedure :: Grid_Update_Host

      ! Procedures to create working arrays to device
      ! (No need to copy them, they are temporrary by nature)
      procedure :: Work_Create_On_Device
      procedure :: Work_Destroy_On_Device

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
#   include "Gpu_Mod/Sparse_Con/Copy_To_Device.f90"
#   include "Gpu_Mod/Sparse_Con/Destroy_On_Device.f90"
#   include "Gpu_Mod/Sparse_Val/Copy_To_Device.f90"
#   include "Gpu_Mod/Sparse_Val/Destroy_On_Device.f90"

    ! Procedures to copy native solver to device
#   include "Gpu_Mod/Native/Copy_To_Device.f90"
#   include "Gpu_Mod/Native/Destroy_On_Device.f90"

    ! Procedures to copy field to device
#   include "Gpu_Mod/Field/Copy_To_Device.f90"
#   include "Gpu_Mod/Field/Destroy_On_Device.f90"
#   include "Gpu_Mod/Field/Update_Host.f90"

    ! Procedures to copy grid to device
#   include "Gpu_Mod/Grid/Copy_To_Device.f90"
#   include "Gpu_Mod/Grid/Destroy_On_Device.f90"
#   include "Gpu_Mod/Grid/Update_Host.f90"

    ! Procedures to create working arrays to device
#   include "Gpu_Mod/Work/Create_On_Device.f90"
#   include "Gpu_Mod/Work/Destroy_On_Device.f90"

  end module
