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
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:)
  real, contiguous, pointer :: v_flux(:), v_m(:), fc(:)
  real                      :: a12, b_tmp, max_abs_val
  real                      :: u_f, v_f, w_f
  real                      :: area_in, area_out, vol_in, vol_out, ratio
  integer                   :: s, c1, c2, i_cel, c, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Insert_Volume_Source_For_Pressure')

  b      => Flow % Nat % b
  v_m    => Flow % Nat % M % v_m
  fc     => Flow % Nat % C % fc
  v_flux => Flow % v_flux

  ! Nullify the volume source
  !$acc kernels
  b(:) = 0.0
  !$acc end kernels

  !------------------------------------------!
  !                                          !
  !   Calculate / update volume flow rates   !
  !                                          !
  !------------------------------------------!

  !----------------------------------------------------!
  !   Calculate volume fluxes through boundary faces   !
  !----------------------------------------------------!
  do reg = Boundary_Regions()

    !$acc parallel loop
    do s = Faces_In_Region(reg)
      c2 = Grid % faces_c(2,s)  ! boundary cell
      v_flux(s) = u_n(c2) * Grid % sx(s)  &
                + v_n(c2) * Grid % sy(s)  &
                + w_n(c2) * Grid % sz(s)
    end do
    !$acc end parallel

  end do

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

      !$acc parallel loop reduction(+:area_in,vol_in)
      do s = Faces_In_Region(reg)
        area_in = area_in + Grid % s(s)
        vol_in  = vol_in  - v_flux(s)
      end do
      !$acc end parallel loop

    end if

    if(Grid % region % type(reg) .eq. OUTFLOW .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$acc parallel loop reduction(+:area_out,vol_out)
      do s = Faces_In_Region(reg)
        area_out = area_out + Grid % s(s)
        vol_out  = vol_out  + v_flux(s)
      end do
      !$acc end parallel

    end if
  end do

  call Global % Sum_Reals(vol_in, area_in)

  !-----------------------------------------------------------------!
  !   If there is volume imbalance in the source for the pressure   !
  !   (correction) equation, pressure will converge poorly, if at   !
  !   all. Thus, correct the outlet fluxes to enforce the balance   !
  !-----------------------------------------------------------------!
  if(.not. Math % Approx_Real(vol_in, vol_out, FEMTO)) then
    ratio = vol_in / vol_out
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW .or.  &
         Grid % region % type(reg) .eq. CONVECT) then

        !$acc parallel loop
        do s = Faces_In_Region(reg)
          v_flux(s) = v_flux(s) * ratio
        end do
        !$acc end parallel

      end if
    end do
  end if

  !--------------------------------------------------!
  !   Calculate volume fluxes through inside faces   !
  !- - - - - - - - - - - - - - - - - - - - - - - - - !
  !   This is application of Rhie & Chow technique   !
  !--------------------------------------------------!

  !$acc parallel loop
  do s = Faces_In_Domain_And_At_Buffers()

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

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
    v_flux(s) = u_f * Grid % sx(s) + v_f * Grid % sy(s) + w_f * Grid % sz(s)  &
              + a12 * (p_n(c1) - p_n(c2))

  end do
  !$acc end parallel

  !---------------------------------------------------------------------!
  !                                                                     !
  !   Calculate volume sources for the pressure (correction) equation   !
  !                                                                     !
  !---------------------------------------------------------------------!

  !---------------------------------!
  !   First consider inside faces   !
  !---------------------------------!

  !$acc parallel loop
  do c1 = Cells_In_Domain()

    b_tmp = b(c1)
    !$acc loop seq
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        b_tmp = b_tmp - v_flux(s) * merge(1,-1, c1.lt.c2)
      end if
    end do
    !$acc end loop

    ! Finish
    b(c1) = b_tmp

  end do
  !$acc end parallel

  !-----------------------------!
  !   Then the boundary faces   !
  !-----------------------------!

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW  .or.  &
       Grid % region % type(reg) .eq. OUTFLOW .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$acc parallel loop
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)  ! inside cell
        b(c1) = b(c1) - v_flux(s)
      end do
      !$acc end parallel

    end if
  end do

  !------------------------------------------------------------------!
  !   Find the cell with the maximum volume imbalance and print it   !
  !------------------------------------------------------------------!
  max_abs_val = 0.0
  !$acc parallel loop reduction(max:max_abs_val)
  do c = Cells_In_Domain()
    max_abs_val = max(max_abs_val, abs(b(c)))
  end do

  ! Find global maximum over all processors
  call Global % Max_Real(max_abs_val)

  !@ Use this for REPORT_VOLUME_BALANCE somehow?
  !@ O_Print '(a,es12.3)', ' # Max. volume balance error '//  &
  !@                       'before correction: ', max_abs_val

  call Profiler % Stop('Insert_Volume_Source_For_Pressure')

  end subroutine
