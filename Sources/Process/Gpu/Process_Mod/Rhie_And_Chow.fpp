!==============================================================================!
  subroutine Rhie_And_Chow(Process, Flow, Grid)
!------------------------------------------------------------------------------!
!>  Second part of inserting volume source for pressure-Poisson equation
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), p_x(:), p_y(:), p_z(:), fc(:)
  real, contiguous, pointer :: u_n(:), v_n(:), w_n(:)
  real                      :: a12, b_tmp
  real                      :: u_f, v_f, w_f, w1, w2
  integer                   :: s, c1, c2, i_cel, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Rhie_And_Chow')

  ! Take some aliases
  ! GPU version doesn't work if you use directly Flow % whatever_variable
  ! These aliases are really needed, not just some gimmick to shorten the code
  b    => Flow % Nat % b
  fc   => Flow % Nat % A % fc

  ! Check if you have pressure gradients at hand and then set aliases properly
  Assert(Flow % stores_gradients_of .eq. 'P')

  ! GPU version doesn't work if you use directly Flow % phi_x, _y and _z
  ! These aliases are really needed, not just some gimmick to shorten the code
  p_x => Flow % phi_x
  p_y => Flow % phi_y
  p_z => Flow % phi_z
  u_n => Flow % u % n
  v_n => Flow % v % n
  w_n => Flow % w % n

  !--------------------------------------------------!
  !   Calculate volume fluxes through inside faces   !
  !--------------------------------------------------!

  !$tf-acc loop begin
  do s = Faces_In_Domain_And_At_Buffers()  ! all present

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    w1 = Grid % f(s)
    w2 = 1.0 - w1

    ! Velocity plus the cell-centered pressure gradient
    ! Units: kg / (m^2 s^2) * m^3 * s / kg = m / s
    u_f = w1 * (u_n(c1) + p_x(c1) * Flow % v_m(c1))  &
        + w2 * (u_n(c2) + p_x(c2) * Flow % v_m(c2))
    v_f = w1 * (v_n(c1) + p_y(c1) * Flow % v_m(c1))  &
        + w2 * (v_n(c2) + p_y(c2) * Flow % v_m(c2))
    w_f = w1 * (w_n(c1) + p_z(c1) * Flow % v_m(c1))  &
        + w2 * (w_n(c2) + p_z(c2) * Flow % v_m(c2))

    ! This is a bit of a code repetition, the
    ! same thing is in the Form_Pressure_Matrix
    ! Anyhow, units are given here are
    !  Units: m * m^3 s / kg = m^4 s / kg
    a12 = fc(s) * (w1 * Flow % v_m(c1) + w2 * Flow % v_m(c2))

    ! Volume flux without the cell-centered pressure gradient
    ! but with the staggered pressure difference
    ! Units:  m^4 s / kg * kg / (m s^2) = m^3 / s
    Flow % v_flux % n(s) = u_f * Grid % sx(s)  &
                         + v_f * Grid % sy(s)  &
                         + w_f * Grid % sz(s)  &
                         + a12 * (Flow % p % n(c1) - Flow % p % n(c2))

  end do
  !$tf-acc loop end

  !---------------------------------------------------------------------!
  !                                                                     !
  !   Calculate volume sources for the pressure (correction) equation   !
  !                                                                     !
  !---------------------------------------------------------------------!

  !---------------------------------!
  !   First consider inside faces   !
  !---------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    b_tmp = b(c1)
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        b_tmp = b_tmp - Flow % v_flux % n(s) * merge(1,-1, c1.lt.c2)
      end if
    end do

    ! Finish
    b(c1) = b_tmp

  end do
  !$tf-acc loop end

  !-----------------------------!
  !   Then the boundary faces   !
  !-----------------------------!

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW  .or.  &
       Grid % region % type(reg) .eq. OUTFLOW .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)   ! inside cell
        b(c1) = b(c1) - Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if
  end do

  call Profiler % Stop('Rhie_And_Chow')

  end subroutine
