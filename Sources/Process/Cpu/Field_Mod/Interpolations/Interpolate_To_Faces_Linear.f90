!==============================================================================!
  subroutine Interpolate_To_Faces_Linear(Flow, phi_f, phi_c,  &
                                         phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Interpolates to all faces in the domain, using the values in cells around  !
!   it, as well as the gradients in these cells.  It is only accurate if the   !
!   cell gradients are computed properly, and the main use of this function    !
!   was intended to be ebmedded in an iterative algorithm for gradient calcu-  !
!   lation by Gaussian theorem.  See also it's sister function Grad_Gauss and  !
!   parent function Grad_Gauss_Variable from this module.                      !
!                                                                              !
!   With OpenMP, this procedure got a speedup of 2.6 on 1M mesh and 4 threads. !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow
  real                      :: phi_f(  Flow % pnt_grid % n_faces)
  real                      :: phi_c( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
  real, optional            :: phi_x( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
  real, optional            :: phi_y( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
  real, optional            :: phi_z( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),     pointer :: Grid
  real,    contiguous, pointer :: xf(:), yf(:), zf(:), xc(:), yc(:), zc(:)
  real,    contiguous, pointer :: rx(:), ry(:), rz(:), f(:)
  integer, contiguous, pointer :: faces_c(:,:)
  integer                      :: s, c1, c2, reg
!==============================================================================!

  call Profiler % Start('Interpolate_To_Faces_Linear')

  ! Take alias
  ! OpenMP doesn't unerstand Fortran's members (%), that's why ...
  ! ... aliases for faces_c, xf, yf, zf, xc ... and rz are needed
  Grid    => Flow % pnt_grid
  faces_c => Grid % faces_c
  f       => Grid % f
  xf      => Grid % xf
  yf      => Grid % yf
  zf      => Grid % zf
  xc      => Grid % xc
  yc      => Grid % yc
  zc      => Grid % zc
  rx      => Grid % rx
  ry      => Grid % ry
  rz      => Grid % rz

  ! Refresh buffers for gradient components (needed)
  call Grid % Exchange_Cells_Real(phi_c)
  call Grid % Exchange_Cells_Real(phi_x)
  call Grid % Exchange_Cells_Real(phi_y)
  call Grid % Exchange_Cells_Real(phi_z)

  !-------------------------------------------------------!
  !   Estimate values at faces from the values in cells   !
  !-------------------------------------------------------!

  ! Perform linear average for boundary faces
  ! (Why doesn't it take care of boundary conditions? - check this!)
  do reg = Boundary_Regions()
    !$omp parallel do     &
    !$omp private(s, c1)  &
    !$omp shared(faces_c, phi_f, phi_c)
    do s = Faces_In_Region(reg)
      c1 = faces_c(1,s)
      phi_f(s) = phi_c(c1)
    end do
    !$omp end parallel do
  end do

  ! Perform linear average for inside faces
  !$omp parallel do                       &
  !$omp private(s, c1, c2)                &
  !$omp shared(faces_c, phi_f, f, phi_c)
  do s = Faces_In_Domain_And_At_Buffers()
    c1  = faces_c(1,s)
    c2  = faces_c(2,s)

    phi_f(s) = f(s) * phi_c(c1) + (1.0 - f(s)) * phi_c(c2)
  end do
  !$omp end parallel do

  !-----------------------------------------------------!
  !   If gradients present, improve the interpolation   !
  !-----------------------------------------------------!
  if(present(phi_x)) then

    do reg = Boundary_Regions()
      !$omp parallel do     &
      !$omp private(s, c1)  &
      !$omp shared(faces_c, phi_f, phi_x, phi_y, phi_z, xf, yf, zf, xc, yc, zc)
      do s = Faces_In_Region(reg)
        c1 = faces_c(1,s)
        phi_f(s) = phi_f(s)                      &
                 + phi_x(c1) * (xf(s) - xc(c1))  &
                 + phi_y(c1) * (yf(s) - yc(c1))  &
                 + phi_z(c1) * (zf(s) - zc(c1))
      end do  ! faces
      !$omp end parallel do
    end do    ! all boundary regions

    !$omp parallel do                                                 &
    !$omp private(s, c1, c2)                                          &
    !$omp shared(faces_c, phi_f, f, phi_x, phi_y, phi_z, rx, ry, rz)
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = faces_c(1,s)
      c2 = faces_c(2,s)

      phi_f(s) = phi_f(s)                       &
               + f(s)                           &
                      * (  phi_x(c1) * rx(s)    &
                         + phi_y(c1) * ry(s)    &
                         + phi_z(c1) * rz(s) )  &
               + (1.0 - f(s))                   &
                      * (  phi_x(c2) * rx(s)    &
                         + phi_y(c2) * ry(s)    &
                         + phi_z(c2) * rz(s) )
    end do
    !$omp end parallel do

  end if  ! gradients are present

  call Profiler % Stop('Interpolate_To_Faces_Linear')

  end subroutine
