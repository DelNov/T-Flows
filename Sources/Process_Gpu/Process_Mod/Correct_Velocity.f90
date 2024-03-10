!==============================================================================!
  subroutine Correct_Velocity(Proc, Flow)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Dimension of the system under consideration                                !
!     [M]{u} = {b}   [kgm/s^2]   [N]                                           !
!                                                                              !
!   Pressure gradient alone:                                                   !
!     p % x        [kg / (m^2 s^2)]                                            !
!                                                                              !
!   Pressure gradient times volume:                                            !
!     p % x * vol  [kg / (m^2 s^2) * m^3 = kg m / s^2 = N]                     !
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  real,    contiguous, pointer :: b(:), v_flux(:), p_n(:), pp_n(:)
  real,    contiguous, pointer :: u_n(:), v_n(:), w_n(:)
  real,    contiguous, pointer :: pp_x(:), pp_y(:), pp_z(:)
  real,    contiguous, pointer :: v_m(:), fc(:)
  integer, contiguous, pointer :: faces_c(:,:)
  integer, contiguous, pointer :: cells_n_cells(:)
  integer, contiguous, pointer :: cells_c(:,:), cells_f(:,:)
  real                         :: a12, b_tmp, max_abs_val
  integer                      :: c, s, c1, c2, nb, nc, nf, i_cel
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Correct_Velocity')

  ! Take some aliases
  b             => Flow % Nat % b
  v_m           => Flow % Nat % M % v_m
  fc            => Flow % Nat % M % fc
  v_flux        => Flow % v_flux
  pp_n          => Flow % pp % n
  p_n           => Flow % p % n
  u_n           => Flow % u % n
  v_n           => Flow % v % n
  w_n           => Flow % w % n
  pp_x          => Flow % pp % x
  pp_y          => Flow % pp % y
  pp_z          => Flow % pp % z
  faces_c       => Flow % pnt_grid % faces_c
  cells_n_cells => Flow % pnt_grid % cells_n_cells
  cells_c       => Flow % pnt_grid % cells_c
  cells_f       => Flow % pnt_grid % cells_f
  nb            =  Flow % pnt_grid % n_bnd_cells
  nc            =  Flow % pnt_grid % n_cells
  nf            =  Flow % pnt_grid % n_faces

  !----------------------!
  !   Correct velocity   !
  !----------------------!

  ! Units: kg m / s^2 * s / kg = m / s
  !$acc parallel loop independent
  do c = 1, nc
    u_n(c) = u_n(c) - pp_x(c) * v_m(c)
    v_n(c) = v_n(c) - pp_y(c) * v_m(c)
    w_n(c) = w_n(c) - pp_z(c) * v_m(c)
  end do
  !$acc end parallel

  !---------------------------------------------!
  !   Correct volume fluxes inside the domain   !
  !---------------------------------------------!

  ! Units: m * m^3 * s / kg * kg / (m s^2) = m^3 / s
  !$acc parallel loop independent
  do s = nb + 1, nf
    c1 = faces_c(1, s)
    c2 = faces_c(2, s)

    a12 = -fc(s) * 0.5 * (v_m(c1) + v_m(c2))

    v_flux(s) = v_flux(s) + (pp_n(c2) - pp_n(c1)) * a12
  end do
  !$acc end parallel

  !-------------------------------!
  !   Re-compute volume sources   !
  !-------------------------------!
  !$acc kernels
  b(:) = 0
  !$acc end kernels

  !$acc parallel loop independent
  do c1 = 1, nc

    b_tmp = b(c1)
    !$acc loop seq
    do i_cel = 1, cells_n_cells(c1)
      c2 = cells_c(i_cel, c1)
      s  = cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        b_tmp = b_tmp - v_flux(s) * merge(1,0, c1.lt.c2)
        b_tmp = b_tmp + v_flux(s) * merge(1,0, c1.gt.c2)
      end if
    end do
    !$acc end loop

    ! Finish, and nullify if it is not in fluid
    b(c1) = b_tmp

  end do
  !$acc end parallel

  !------------------------------------------------------------------!
  !   Find the cell with the maximum volume imbalance and print it   !
  !------------------------------------------------------------------!
  max_abs_val = 0.0
  !$acc parallel loop reduction(max:max_abs_val)
  do c = 1, nc
    max_abs_val = max(max_abs_val, abs(b(c)))
  end do
  print '(a,es12.3)', ' # Max. volume balance error '//  &
                      'after correction: ', max_abs_val

  !-----------------------------------!
  !     Update the pressure field     !
  !   (hard-coded under-relaxation)   !
  !-----------------------------------!

  !$acc parallel loop independent
  do c = 1, nc
    p_n(c) = p_n(c) + 0.2 * pp_n(c)
  end do
  !$acc end parallel

  call Profiler % Stop('Correct_Velocity')

  end subroutine
