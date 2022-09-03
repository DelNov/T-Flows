!==============================================================================!
!   Introduce new types to be used with User_Mod                               !
!==============================================================================!

  !-----------------!
  !                 !
  !   Tensor type   !
  !                 !
  !-----------------!
  type Tensor_Type

    real, allocatable :: x (:)
    real, allocatable :: y (:)
    real, allocatable :: z (:)
    real, allocatable :: xy(:)
    real, allocatable :: xz(:)
    real, allocatable :: yz(:)

  end type

  type(Tensor_Type) :: cell_inertia
