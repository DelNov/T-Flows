!==============================================================================!
  subroutine Insert_Volume_Source_For_Pressure(Proc, Flow)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Discreetized system of momentum conservation equations:                    !
!                                                                              !
!     [M]{u} = {b}   [kg m/s^2]   [N]                                          !
!                                                                              !
!   Dimensions of certain variables:                                           !
!                                                                              !
!     M               [kg/s]                                                   !
!     u, v, w         [m/s]                                                    !
!     b               [kg m/s^2]     [N]                                       !
!     p, pp           [kg/(m s^2)]   [N/m^2]                                   !
!     p%x, p%y, p%z   [kg/(m^2 s^2)]                                           !
!     v_flux          [m^3/s]                                                  !
!------------------------------------------------------------------------------!
!   Discretized pressure-Poisson equation reads:                               !
!                                                                              !
!     [A] {pp} = {b}     [m^3/s]                                               !
!                                                                              !
!   Dimensions of certain variables:                                           !
!                                                                              !
!     A               [m^4 s/kg]                                               !
!     pp              [kg/(m s^2)]                                             !
!     p%x, p%y, p%z   [kg/(m^2 s^2)]                                           !
!     b               [m^3/s]                                                  !
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  real,    contiguous, pointer :: b(:), u_n(:), v_n(:), w_n(:)
  real,    contiguous, pointer :: p_n(:), p_x(:), p_y(:), p_z(:)
  real,    contiguous, pointer :: v_flux(:), v_m(:), sx(:), sy(:), sz(:), fc(:)
  integer, contiguous, pointer :: faces_c(:,:)
  integer, contiguous, pointer :: cells_n_cells(:)
  integer, contiguous, pointer :: cells_c(:,:), cells_f(:,:)
  real                         :: a12, b_tmp, max_abs_val
  real                         :: u_f, v_f, w_f
  integer                      :: s, c1, c2, nb, nc, nf, i_cel, c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Insert_Volume_Source_For_Pressure')

  b             => Flow % Nat % b
  v_m           => Flow % Nat % M % v_m
  fc            => Flow % Nat % M % fc
  p_n           => Flow % p % n
  p_x           => Flow % p % x
  p_y           => Flow % p % y
  p_z           => Flow % p % z
  u_n           => Flow % u % n
  v_n           => Flow % v % n
  w_n           => Flow % w % n
  v_flux        => Flow % v_flux
  faces_c       => Flow % pnt_grid % faces_c
  nb            =  Flow % pnt_grid % n_bnd_cells
  nc            =  Flow % pnt_grid % n_cells
  nf            =  Flow % pnt_grid % n_faces
  sx            => Flow % pnt_grid % sx
  sy            => Flow % pnt_grid % sy
  sz            => Flow % pnt_grid % sz
  cells_n_cells => Flow % pnt_grid % cells_n_cells
  cells_c       => Flow % pnt_grid % cells_c
  cells_f       => Flow % pnt_grid % cells_f

  ! Nullify the volume source
  !$acc kernels
  b(:) = 0.0
  !$acc end kernels

  !-------------------------------------------!
  !   Calculate volume fluxes through faces   !
  !-------------------------------------------!

  !$acc parallel loop
  do s = nb + 1, nf

    c1 = faces_c(1,s)
    c2 = faces_c(2,s)

    ! Velocity plus the cell-centered pressure gradient
    ! Units: kg / (m^2 s^2) * m^3 * s / kg = m / s
    u_f = 0.5 * (u_n(c1) + p_x(c1) * v_m(c1) + u_n(c2) + p_x(c2) * v_m(c2))
    v_f = 0.5 * (v_n(c1) + p_y(c1) * v_m(c1) + v_n(c2) + p_y(c2) * v_m(c2))
    w_f = 0.5 * (w_n(c1) + p_z(c1) * v_m(c1) + w_n(c2) + p_z(c2) * v_m(c2))

    ! This is a bit of a code repetition, the
    ! same thing is in the Form_Pressure_Matrix
    ! Anyhow, units are given here are
    !  Units: m * m^3 s / kg = m^4 s / kg
    a12 = fc(s) * 0.5 * (v_m(c1) + v_m(c2))

    ! Volume flux without the cell-centered pressure gradient
    ! but with the staggered pressure difference
    ! Units:  m^4 s / kg * kg / (m s^2) = m^3 / s
    v_flux(s) = u_f * sx(s) + v_f * sy(s) + w_f * sz(s)  &
              + a12 * (p_n(c1) - p_n(c2))

  end do
  !$acc end parallel

  !----------------------------------------------------!
  !   Calculate volume sources with corrected fluxes   !
  !----------------------------------------------------!

  !$acc parallel loop
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

    ! Finish
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
                      'before correction: ', max_abs_val

  call Profiler % Stop('Insert_Volume_Source_For_Pressure')

  end subroutine
