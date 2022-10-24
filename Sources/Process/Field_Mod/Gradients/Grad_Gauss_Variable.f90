!==============================================================================!
  subroutine Grad_Gauss_Variable(Flow, phi)
!------------------------------------------------------------------------------!
!   Tries to find gradients with Gaussian in an iterative fashion.  It works   !
!   really bad for tetrahedral grids with no initial values of gradients, and  !
!   a more elaborate approach is therefore needed, which will probably be in   !
!   Grad_Gauss_Pressure, when introduced.
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow
  type(Var_Type),    target :: phi
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer  :: Grid
  integer                   :: s, c, c1, c2, iter
  real                      :: res, norm
  real, contiguous, pointer :: phi_f_n(:), phi_f_o(:)
!==============================================================================!

  call Profiler % Start('Grad_Gauss_Variable')

  call Work % Connect_Real_Face(phi_f_n, phi_f_o)

  ! Take alias
  Grid => Flow % pnt_grid

  phi_f_n(:) = 0.0
  phi_f_o(:) = 0.0

  norm = max(maxval(phi % n) - minval(phi % n), MICRO)
  call Comm_Mod_Global_Max_Real(norm)

  !-------------------------------!
  !   Start iterative procedure   !
  !-------------------------------!
  do iter = 1, Flow % gauss_miter

    ! Save the old iteration
    phi_f_o(:) = phi_f_n(:)

    ! Estimate values at faces from the
    ! values in cells and last gradients
    call Flow % Interpolate_To_Faces_Linear(phi_f_n, phi % n,  &
                                            phi % x, phi % y, phi % z)

    ! Update gradients from the values at faces
    call Flow % Grad_Gauss(phi_f_n, phi % x, phi % y, phi % z)

    ! Take the difference between two iterations
    res = maxval(abs(phi_f_n(:)-phi_f_o(:))) / norm
    call Comm_Mod_Global_Max_Real(res)

    if(res < Flow % gauss_tol) exit

  end do  ! iter

  ! Accumulate all Gauss iterations
  Flow % gauss_iters = Flow % gauss_iters + iter
  Flow % gauss_calls = Flow % gauss_calls + 1

  ! Refresh buffers for gradient components (these calls are needed)
  call Grid % Exchange_Cells_Real(phi % x)
  call Grid % Exchange_Cells_Real(phi % y)
  call Grid % Exchange_Cells_Real(phi % z)

  ! if(this_proc < 2) then
  !   print *, '# Final residual from Gauss: ', res,  &
  !            ' reached in ', iter, ' iterations '
  ! end if

  call Work % Disconnect_Real_Face(phi_f_n, phi_f_o)

  call Profiler % Stop('Grad_Gauss_Variable')

  end subroutine
