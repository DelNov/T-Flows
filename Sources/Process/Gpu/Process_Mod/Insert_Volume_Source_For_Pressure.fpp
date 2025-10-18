!==============================================================================!
  subroutine Insert_Volume_Source_For_Pressure(Process, Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Discreetized system of momentum conservation equations:                    !
!                                                                              !
!     [M]{u} = {b}   [kg m/s^2]   [N]                                          !
!                                                                              !
!   Dimensions of certain variables:                                           !
!                                                                              !
!     M                    [kg/s]                                              !
!     u, v, w              [m/s]                                               !
!     b                    [kg m/s^2]     [N]                                  !
!     p, pp                [kg/(m s^2)]   [N/m^2]                              !
!     p % x, p % y, p % z  [kg/(m^2 s^2)]                                      !
!     v_flux % n           [m^3/s]                                             !
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
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), p_x(:), p_y(:), p_z(:), fc(:)
  real, contiguous, pointer :: u_n(:), v_n(:), w_n(:)
  real                      :: vol_in, vol_out, area_in, area_out
  real                      :: a12, b_tmp, fac, max_abs_val
  real                      :: u_f, v_f, w_f, w1, w2
  integer                   :: s, c1, c2, i_cel, c, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Insert_Volume_Source_For_Pressure')

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

  ! Nullify the volume source
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()  ! all present
    b(c) = 0.0
  end do
  !$tf-acc loop end

  !------------------------------------------!
  !                                          !
  !   Calculate / update volume flow rates   !
  !                                          !
  !------------------------------------------!

  !----------------------------------------------------!
  !   Calculate volume fluxes through boundary faces   !
  !----------------------------------------------------!
  do reg = Boundary_Regions()

    if(Grid % region % type(reg) .eq. SYMMETRY) then
      !$tf-acc loop begin
      do s = Faces_In_Region(reg)
        Flow % v_flux % n(s) = 0.0
      end do
      !$tf-acc loop end
    else
      !$tf-acc loop begin
      do s = Faces_In_Region(reg)
        c2 = Grid % faces_c(2,s)  ! boundary cell
        Flow % v_flux % n(s) = Flow % u % n(c2) * Grid % sx(s)  &
                             + Flow % v % n(c2) * Grid % sy(s)  &
                             + Flow % w % n(c2) * Grid % sz(s)
      end do  ! faces
      !$tf-acc loop end
    end if    ! boundary region type

  end do      ! regions

  !-------------------------------------------!
  !   Calcululate inflow and outflow volume   !
  !   flow rates and inlet and outlet areas   !
  !-------------------------------------------!
  vol_in   = 0.0
  vol_out  = 0.0
  area_in  = 0.0
  area_out = 0.0

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        area_in = area_in + Grid % s(s)
        vol_in  = vol_in  - Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if

    if(Grid % region % type(reg) .eq. PRESSURE .or.  &
       Grid % region % type(reg) .eq. OUTFLOW  .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        area_out = area_out + Grid % s(s)
        vol_out  = vol_out  + Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if
  end do

  call Global % Sum_Reals(vol_in,   &
                          vol_out,  &
                          area_in,  &
                          area_out)

  !-----------------------------------------------------------------!
  !   If there is volume imbalance in the source for the pressure   !
  !   (correction) equation, pressure will converge poorly, if at   !
  !   all. Thus, correct the outlet fluxes to enforce the balance   !
  !-----------------------------------------------------------------!

  if(.not. Math % Approx_Real(vol_in, vol_out, FEMTO)) then

    ! You should compute the "fac" now ...
    fac = vol_in / (vol_out + TINY)

    ! ... and correct all velocities
    vol_out = 0.0

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW   .or.  &
         Grid % region % type(reg) .eq. CONVECT   .or.  &
         Grid % region % type(reg) .eq. PRESSURE) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)  ! all present
          c2 = Grid % faces_c(2,s)

          ! Update velocity components ...
          Flow % u % n(c2) = Flow % u % n(c2) * fac
          Flow % v % n(c2) = Flow % v % n(c2) * fac
          Flow % w % n(c2) = Flow % w % n(c2) * fac

          ! ... volume flux itself ...
          Flow % v_flux % n(s) = Flow % u % n(c2) * Grid % sx(s)    &
                               + Flow % v % n(c2) * Grid % sy(s)    &
                               + Flow % w % n(c2) * Grid % sz(s)

          ! ... and bulk volume out
          vol_out = vol_out + Flow % v_flux % n(s)
        end do
        !$tf-acc loop end

      end if
    end do
  end if

  ! Holy mackrele: summ it up over all processors
  call Global % Sum_Real(vol_out)  ! not checked

  !--------------------------------------------------!
  !   Calculate volume fluxes through inside faces   !
  !- - - - - - - - - - - - - - - - - - - - - - - - - !
  !   This is application of Rhie & Chow technique   !
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

  !------------------------------------------------------------------!
  !   Find the cell with the maximum volume imbalance and print it   !
  !------------------------------------------------------------------!
  max_abs_val = 0.0
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    max_abs_val = max(max_abs_val, abs(b(c)))
  end do
  !$tf-acc loop end

  ! Find global maximum over all processors
  call Global % Max_Real(max_abs_val)

  !@ Use this for REPORT_VOLUME_BALANCE somehow?
  !@ O_Print '(a,es12.3)', ' # Max. volume balance error '//  &
  !@                       'before correction: ', max_abs_val

  call Profiler % Stop('Insert_Volume_Source_For_Pressure')

  end subroutine
