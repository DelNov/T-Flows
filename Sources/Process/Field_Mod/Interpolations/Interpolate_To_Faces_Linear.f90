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
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c1, c2, reg
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

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
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      phi_f(s) = phi_c(c1)
    end do
  end do

  ! Perform linear average for inside faces
  do s = Faces_In_Domain()
    c1  = Grid % faces_c(1,s)
    c2  = Grid % faces_c(2,s)

    phi_f(s) = Grid % f(s) * phi_c(c1) + (1.0 - Grid % f(s)) * phi_c(c2)
  end do

  !-----------------------------------------------------!
  !   If gradients present, improve the interpolation   !
  !-----------------------------------------------------!
  if(present(phi_x)) then

    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        phi_f(s) = phi_f(s)                                    &
                 + phi_x(c1) * (Grid % xf(s) - Grid % xc(c1))  &
                 + phi_y(c1) * (Grid % yf(s) - Grid % yc(c1))  &
                 + phi_z(c1) * (Grid % zf(s) - Grid % zc(c1))
      end do  ! faces
    end do    ! all boundary regions

    do s = Faces_In_Domain()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      phi_f(s) = phi_f(s)                                     &
               + Grid % f(s)                                  &
                             * (  phi_x(c1) * Grid % rx(s)    &
                                + phi_y(c1) * Grid % ry(s)    &
                                + phi_z(c1) * Grid % rz(s) )  &
               + (1.0 - Grid % f(s))                          &
                             * (  phi_x(c2) * Grid % rx(s)    &
                                + phi_y(c2) * Grid % ry(s)    &
                                + phi_z(c2) * Grid % rz(s) )
    end do

  end if  ! gradients are present

  end subroutine
