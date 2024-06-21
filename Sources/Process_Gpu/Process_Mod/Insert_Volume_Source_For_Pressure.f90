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
  type(Grid_Type)          :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), p_x(:), p_y(:), p_z(:), fc(:)
  real, contiguous, pointer :: u_n(:), v_n(:), w_n(:)
  real                      :: a12, b_tmp, max_abs_val
  real                      :: u_f, v_f, w_f
  real                      :: area_in, area_out, vol_in, vol_out, ratio
  integer                   :: s, c1, c2, i_cel, c, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Insert_Volume_Source_For_Pressure')

  ! Take some aliases
  ! GPU version doesn't work if you use directly Flow % whatever_variable
  ! These aliases are really needed, not just some gimmick to shorten the code
  b        => Flow % Nat % b
  fc       => Flow % Nat % C % fc

  ! Check if you have pressure gradients at hand and then set aliases properly
  Assert(Flow % stores_gradients_of .eq. 'P')

  ! GPU version doesn't work if you use directly Flow % phi_x, _y and _z
  ! These aliases are really needed, not just some gimmick to shorten the code
  p_x => Flow % phi_x
  p_y => Flow % phi_y
  p_z => Flow % phi_z
  u_n => flow_u_n
  v_n => flow_v_n
  w_n => flow_w_n

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

    !$acc parallel loop                                                  &
    !$acc present(grid_faces_c, grid_region_f_face, grid_region_l_face,  &
    !$acc         grid_sx, grid_sy, grid_sz,                             &
    !$acc         u_n, v_n, w_n, flow_v_flux_n)
    do s = Faces_In_Region_Gpu(reg)  ! all present
      c2 = grid_faces_c(2,s)  ! boundary cell
      flow_v_flux_n(s) = u_n(c2) * grid_sx(s)  &
                       + v_n(c2) * grid_sy(s)  &
                       + w_n(c2) * grid_sz(s)
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

      !$acc parallel loop reduction(+:area_in,vol_in)        &
      !$acc present(grid_region_f_face, grid_region_l_face,  &
      !$acc         grid_s,                                  &
      !$acc         flow_v_flux_n)
      do s = Faces_In_Region_Gpu(reg)  ! all present
        area_in = area_in + grid_s(s)
        vol_in  = vol_in  - flow_v_flux_n(s)
      end do
      !$acc end parallel loop

    end if

    if(Grid % region % type(reg) .eq. OUTFLOW .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$acc parallel loop reduction(+:area_out,vol_out)      &
      !$acc present(grid_region_f_face, grid_region_l_face,  &
      !$acc         grid_s,                                  &
      !$acc         flow_v_flux_n)
      do s = Faces_In_Region_Gpu(reg)  ! all present
        area_out = area_out + grid_s(s)
        vol_out  = vol_out  + flow_v_flux_n(s)
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

        !$acc parallel loop                                    &
        !$acc present(grid_region_f_face, grid_region_l_face,  &
        !$acc         flow_v_flux_n)
        do s = Faces_In_Region_Gpu(reg)  ! all present
          flow_v_flux_n(s) = flow_v_flux_n(s) * ratio
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

  !$acc parallel loop                                     &
  !$acc present(grid_faces_c,                             &
  !$acc         grid_region_f_face, grid_region_l_face,   &
  !$acc         grid_sx, grid_sy, grid_sz,                &
  !$acc         u_n, v_n, w_n, p_x, p_y, p_z,             &
  !$acc         flow_v_m, flow_v_flux_n, flow_p_n)
  do s = Faces_In_Domain_And_At_Buffers_Gpu()  ! all present

    c1 = grid_faces_c(1,s)
    c2 = grid_faces_c(2,s)

    ! Velocity plus the cell-centered pressure gradient
    ! Units: kg / (m^2 s^2) * m^3 * s / kg = m / s
    u_f = Face_Value(s, u_n(c1) + p_x(c1)*flow_v_m(c1),  u_n(c2) + p_x(c2)*flow_v_m(c2))
    v_f = Face_Value(s, v_n(c1) + p_y(c1)*flow_v_m(c1),  v_n(c2) + p_y(c2)*flow_v_m(c2))
    w_f = Face_Value(s, w_n(c1) + p_z(c1)*flow_v_m(c1),  w_n(c2) + p_z(c2)*flow_v_m(c2))


    ! This is a bit of a code repetition, the
    ! same thing is in the Form_Pressure_Matrix
    ! Anyhow, units are given here are
    !  Units: m * m^3 s / kg = m^4 s / kg
    a12 = fc(s) * Face_Value(s, flow_v_m(c1), flow_v_m(c2))

    ! Volume flux without the cell-centered pressure gradient
    ! but with the staggered pressure difference
    ! Units:  m^4 s / kg * kg / (m s^2) = m^3 / s
    flow_v_flux_n(s) = u_f * grid_sx(s)  &
                     + v_f * grid_sy(s)  &
                     + w_f * grid_sz(s)  &
                     + a12 * (flow_p_n(c1) - flow_p_n(c2))

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

  !$acc parallel loop                                            &
  !$acc present(grid_cells_n_cells, grid_cells_c, grid_cells_f,  &
  !$acc         grid_region_f_cell, grid_region_l_cell,          &
  !$acc         flow_v_flux_n, b)
  do c1 = Cells_In_Domain_Gpu()  ! all present

    b_tmp = b(c1)
    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        b_tmp = b_tmp - flow_v_flux_n(s) * merge(1,-1, c1.lt.c2)
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

      !$acc parallel loop                                                  &
      !$acc present(grid_faces_c, grid_region_f_face, grid_region_l_face,  &
      !$acc         flow_v_flux_n, b)
      do s = Faces_In_Region_Gpu(reg)  ! all present
        c1 = grid_faces_c(1,s)  ! inside cell
        b(c1) = b(c1) - flow_v_flux_n(s)
      end do
      !$acc end parallel

    end if
  end do

  !------------------------------------------------------------------!
  !   Find the cell with the maximum volume imbalance and print it   !
  !------------------------------------------------------------------!
  max_abs_val = 0.0
  !$acc parallel loop reduction(max:max_abs_val)         &
  !$acc present(grid_region_f_cell, grid_region_l_cell,  &
  !$acc         b)
  do c = Cells_In_Domain_Gpu()  ! all present
    max_abs_val = max(max_abs_val, abs(b(c)))
  end do

  ! Find global maximum over all processors
  call Global % Max_Real(max_abs_val)

  !@ Use this for REPORT_VOLUME_BALANCE somehow?
  !@ O_Print '(a,es12.3)', ' # Max. volume balance error '//  &
  !@                       'before correction: ', max_abs_val

  call Profiler % Stop('Insert_Volume_Source_For_Pressure')

  end subroutine
