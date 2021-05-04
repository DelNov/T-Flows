!==============================================================================!
  subroutine Field_Mod_Interpolate_To_Faces_Linear(flow, phi_f, phi_c,  &
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
  type(Field_Type), target :: flow
  real                     :: phi_f(  flow % pnt_grid % n_faces)
  real                     :: phi_c( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
  real, optional           :: phi_x( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
  real, optional           :: phi_y( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
  real, optional           :: phi_z( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, phi_c)
  call Grid_Mod_Exchange_Cells_Real(grid, phi_x)
  call Grid_Mod_Exchange_Cells_Real(grid, phi_y)
  call Grid_Mod_Exchange_Cells_Real(grid, phi_z)

  !-------------------------------------------------------!
  !   Estimate values at faces from the values in cells   !
  !------------------------------------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Inside
    if(c2 > 0) then
      phi_f(s) = grid % f(s)         * phi_c(c1)  &
               + (1.0 - grid % f(s)) * phi_c(c2)

    ! On the boundary cells
    else
      phi_f(s) = phi_c(c1)

    end if

  end do

  !-----------------------------------------------------!
  !   If gradients present, improve the interpolation   !
  !-----------------------------------------------------!
  if(present(phi_x)) then

    do s = 1, grid % n_faces

      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Inside
      if(c2 > 0) then
        phi_f(s) = phi_f(s)                                     &
                 + grid % f(s)                                  &
                               * (  phi_x(c1) * grid % rx(s)    &
                                  + phi_y(c1) * grid % ry(s)    &
                                  + phi_z(c1) * grid % rz(s) )  &
                 + (1.0 - grid % f(s))                          &
                               * (  phi_x(c2) * grid % rx(s)    &
                                  + phi_y(c2) * grid % ry(s)    &
                                  + phi_z(c2) * grid % rz(s) )
      ! On the boundary cells
      else
        phi_f(s) = phi_f(s)                                    &
                 + phi_x(c1) * (grid % xf(s) - grid % xc(c1))  &
                 + phi_y(c1) * (grid % yf(s) - grid % yc(c1))  &
                 + phi_z(c1) * (grid % zf(s) - grid % zc(c1))
      end if

    end do

  end if

  end subroutine
