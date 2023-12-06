!==============================================================================!
  subroutine Grad_Gauss_Variable(Flow, phi)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to calculate gradients of a generic variable
!>  (denoted as phi) using Gauss's theorem. It is an iterative procedure that
!>  continuously refines the gradient calculation for improved accuracy.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Gradient calculation: The subroutine calculates the gradients of a       !
!     specified variable within the flow field using Gauss's theorem.          !
!   * Iterative refinement: The process involves an iterative refinement       !
!     approach to ensure accuracy, particularly in complex flow fields.        !
!   * Applicability: While primarily designed for general variables, it        !
!     is still used from withing the Grad_Gauss_Pressure inside an iterative   !
!     process of computing the gradients of pressure and pressure corrections. !
!------------------------------------------------------------------------------!
!   Notes                                                                      !
!                                                                              !
!   * It works really badly for tetrahedral grids with no initial values of    !
!     gradients, and  a more elaborate approach is therefore needed, probably  !
!     similar to the one already implemented in Grad_Gauss_Pressure.           !
!   * With OpenMP, this procedure got a speedup of 1.6 on 1M mesh & 4 threads. !
!     (I reckon there is not so much potential for improving this, few loops.  !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow object
  type(Var_Type),    target :: phi   !! variable whose gradients are calculated
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer  :: Grid
  integer                   :: s, iter
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
  call Global % Max_Real(norm)

  !-------------------------------!
  !   Start iterative procedure   !
  !-------------------------------!
  do iter = 1, Flow % gauss_miter

    ! Estimate values at faces from the
    ! values in cells and last gradients
    call Flow % Interpolate_To_Faces_Linear(phi_f_n, phi % n,  &
                                            phi % x, phi % y, phi % z)

    ! Update gradients from the values at faces
    call Flow % Grad_Gauss(Grid, phi_f_n, phi % x, phi % y, phi % z)

    ! Take the difference between two iterations to find residual
    ! and copy the new value to the old value along the way
    res = 0.0
    !$omp parallel do private(s) shared(phi_f_n, phi_f_o) reduction(max : res)
    do s = 1, Grid % n_faces
      res = max(res, abs(phi_f_n(s) - phi_f_o(s)))
      phi_f_o(s) = phi_f_n(s)
    end do
    !$omp end parallel do

    res = res / norm
    call Global % Max_Real(res)

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
