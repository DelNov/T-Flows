!==============================================================================!
  subroutine Field_Mod_Grad_Gauss_Variable(flow, phi)
!------------------------------------------------------------------------------!
!   Tries to find gradients with Gaussian in an iterative fashion.  It works   !
!   really bad for tetrahedral grids with no initial values of gradients, and  !
!   a more elaborate approach is therefore needed, which will probably be in   !
!   Grad_Gauss_Pressure, when introduced.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_f_n  => r_face_01,  &
                      phi_f_o  => r_face_02
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type), target :: flow
  type(Var_Type),   target :: phi
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2, iter
  real                     :: res, norm
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  phi_f_n(:) = 0.0
  phi_f_o(:) = 0.0

  norm = max(maxval(phi % n) - minval(phi % n), MICRO)
  call Comm_Mod_Global_Max_Real(norm)

  !-------------------------------!
  !   Start iterative procedure   !
  !-------------------------------!
  do iter = 1, flow % gauss_miter

    ! Save the old iteration
    phi_f_o(:) = phi_f_n(:)

    ! Estimate values at faces from the
    ! values in cells and last gradients
    call Field_Mod_Interpolate_To_Faces(flow, phi_f_n, phi % n,  &
                                              phi % x, phi % y, phi % z)

    ! Update gradients from the values at faces
    call Field_Mod_Grad_Gauss(flow, phi_f_n, phi % x, phi % y, phi % z)

    ! Take the difference between two iterations
    res = maxval(abs(phi_f_n(:)-phi_f_o(:))) / norm
    call Comm_Mod_Global_Max_Real(res)

    if(res < flow % gauss_tol) exit

  end do  ! iter

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, phi % x)
  call Grid_Mod_Exchange_Cells_Real(grid, phi % y)
  call Grid_Mod_Exchange_Cells_Real(grid, phi % z)

  if(this_proc < 2) then
    print *, '# Final residual from Gauss: ', res,  &
             ' reached in ', iter, ' iterations '
  end if

  end subroutine
