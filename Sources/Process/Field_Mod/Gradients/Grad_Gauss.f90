!==============================================================================!
  subroutine Grad_Gauss(Flow, Grid, phi_f, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!>  Calculates gradient of generic variable phi by the Gauss' theorem
!>  and stores its three components in arrays phi_x, phi_y and phi_z.
!>  This subroutine refreshes the buffers of the variable before calculating
!>  the gradients, and of the calculated gradient components.
!------------------------------------------------------------------------------!
!   Notes                                                                      !
!                                                                              !
!   * It heavily relies on the accuracy of values in faces, which are not      !
!     calculated here, and the primary use of this function is to be used as   !
!     embedded in an iterative algorithm which also updates face values.       !
!   * See also it's sister function Interpolate_To_Faces from this module,     !
!     and its parent function Grad_Gauss_Variable, also from this module.      !
!   * With OpenMP, this procedure got a speedup of 3.5 on 1M mesh & 4 threads. !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow object
  type(Grid_Type),   target :: Grid  !! grid object
  real                      :: phi_f(1:Grid % n_faces)
    !! face values of the field whose gradients are being calculated
  real                      :: phi_x( -Grid % n_bnd_cells:Grid % n_cells)
    !! x component of the calculated gradient
  real                      :: phi_y( -Grid % n_bnd_cells:Grid % n_cells)
    !! y component of the calculated gradient
  real                      :: phi_z( -Grid % n_bnd_cells:Grid % n_cells)
    !! z component of the calculated gradient
!----------------------------------[Locals]------------------------------------!
  real,    contiguous, pointer :: sx(:), sy(:), sz(:), vol(:)
  integer, contiguous, pointer :: faces_c(:,:)
  integer                      :: s, c1, c2, c, reg
!==============================================================================!

  call Profiler % Start('Grad_Gauss')

  ! Take alias
  ! OpenMP doesn't unerstand Fortran's members (%), that's ...
  ! ... why aliases for faces_c, sx, sy, sz and vol are needed
  faces_c => Grid % faces_c
  sx      => Grid % sx
  sy      => Grid % sy
  sz      => Grid % sz
  vol     => Grid % vol

  !-----------------------------------------------!
  !   Update gradients from the values at faces   !
  !-----------------------------------------------!
  !$omp parallel do private(c) shared(phi_x, phi_y, phi_z)
  do c = Cells_In_Domain()
    phi_x(c) = 0.0
    phi_y(c) = 0.0
    phi_z(c) = 0.0
  end do
  !$omp end parallel do

  ! Faces in boundary region first
  do reg = Boundary_Regions()
    !$omp parallel do                                              &
    !$omp private(s, c1)                                           &
    !$omp shared(faces_c, phi_f, phi_x, phi_y, phi_z, sx, sy, sz)
    do s = Faces_In_Region(reg)
      c1 = faces_c(1,s)

      phi_x(c1) = phi_x(c1) - phi_f(s) * sx(s)
      phi_y(c1) = phi_y(c1) - phi_f(s) * sy(s)
      phi_z(c1) = phi_z(c1) - phi_f(s) * sz(s)
    end do
    !$omp end parallel do
  end do

  ! Faces inside the domain
  !$omp parallel do                                              &
  !$omp private(s, c1, c2)                                       &
  !$omp shared(faces_c, phi_f, phi_x, phi_y, phi_z, sx, sy, sz)
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = faces_c(1,s)
    c2 = faces_c(2,s)

    phi_x(c1) = phi_x(c1) - phi_f(s) * sx(s)
    phi_y(c1) = phi_y(c1) - phi_f(s) * sy(s)
    phi_z(c1) = phi_z(c1) - phi_f(s) * sz(s)

    phi_x(c2) = phi_x(c2) + phi_f(s) * sx(s)
    phi_y(c2) = phi_y(c2) + phi_f(s) * sy(s)
    phi_z(c2) = phi_z(c2) + phi_f(s) * sz(s)
  end do
  !$omp end parallel do

  !$omp parallel do private(c) shared(phi_x, phi_y, phi_z, vol)
  do c = Cells_In_Domain()
    phi_x(c) = -phi_x(c) / vol(c)
    phi_y(c) = -phi_y(c) / vol(c)
    phi_z(c) = -phi_z(c) / vol(c)
  end do
  !$omp end parallel do

  call Profiler % Stop('Grad_Gauss')

  end subroutine
