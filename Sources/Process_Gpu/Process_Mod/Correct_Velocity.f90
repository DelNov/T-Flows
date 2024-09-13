!==============================================================================!
  subroutine Correct_Velocity(Process, Grid, Flow)
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
  class(Process_Type)      :: Process
  type(Grid_Type)          :: Grid
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), fc(:), pp_x(:), pp_y(:), pp_z(:)
  real, contiguous, pointer :: visc(:), dens(:)
  real                      :: a12, b_tmp, max_abs_val
  real                      :: cfl_max, pe_max, cfl_t, pe_t, nu_f
  integer                   :: c, s, c1, c2, i_cel, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Correct_Velocity')

  ! Take some aliases
  ! GPU version doesn't work if you use directly Flow % whatever_variable
  ! These aliases are really needed, not just some gimmick to shorten the code
  b    => Flow % Nat % b
  fc   => Flow % Nat % C % fc
  dens => flow_density
  visc => flow_viscosity

  ! Check if you have pressure gradients at hand and then set aliases properly
  Assert(Flow % stores_gradients_of .eq. 'PP')

  ! GPU version doesn't work if you use directly Flow % phi_x, _y and _z
  ! These aliases are really needed, not just some gimmick to shorten the code
  pp_x => Flow % phi_x
  pp_y => Flow % phi_y
  pp_z => Flow % phi_z

  !----------------------!
  !                      !
  !   Correct velocity   !
  !                      !
  !----------------------!

  ! Units: kg m / s^2 * s / kg = m / s
  !$acc parallel loop  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   flow_u_n,  &
  !$acc   pp_x,  &
  !$acc   flow_v_m,  &
  !$acc   flow_v_n,  &
  !$acc   pp_y,  &
  !$acc   flow_w_n,  &
  !$acc   pp_z   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present
    flow_u_n(c) = flow_u_n(c) - pp_x(c) * flow_v_m(c)
    flow_v_n(c) = flow_v_n(c) - pp_y(c) * flow_v_m(c)
    flow_w_n(c) = flow_w_n(c) - pp_z(c) * flow_v_m(c)
  end do
  !$acc end parallel

  ! Update buffers for velocities over all processors
  call Grid % Exchange_Cells_Real(Flow % u % n)
  call Grid % Exchange_Cells_Real(Flow % v % n)
  call Grid % Exchange_Cells_Real(Flow % w % n)

  !---------------------------------------------!
  !                                             !
  !   Correct volume fluxes inside the domain   !
  !                                             !
  !---------------------------------------------!

  ! Units: m * m^3 * s / kg * kg / (m s^2) = m^3 / s
  !$acc parallel loop  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   grid_faces_c,  &
  !$acc   fc,  &
  !$acc   flow_v_m,  &
  !$acc   flow_v_flux_n,  &
  !$acc   flow_pp_n   &
  !$acc )
  do s = grid_region_f_face(grid_n_regions), grid_region_l_face(grid_n_regions)  ! all present
    c1 = grid_faces_c(1, s)
    c2 = grid_faces_c(2, s)

    a12 = -fc(s) * Face_Value(s, flow_v_m(c1), flow_v_m(c2))

    flow_v_flux_n(s) = flow_v_flux_n(s)  &
                         + (flow_pp_n(c2) - flow_pp_n(c1)) * a12
  end do
  !$acc end parallel

  !-------------------------------!
  !                               !
  !   Re-compute volume sources   !
  !                               !
  !- - - - - - - - - - - - - - - -!
  !   This step is for checking   !
  !-------------------------------!
  !$acc kernels
  b(:) = 0
  !$acc end kernels

  !---------------------------------!
  !   First consider inside faces   !
  !---------------------------------!

  !$acc parallel loop independent                                &
  !$acc present(grid_cells_c, grid_cells_f,                      &
  !$acc         grid_cells_n_cells, grid_cells_c, grid_cells_f,  &
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

    ! Finish, and nullify if it is not in fluid
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

      !$acc parallel loop  &
      !$acc present(  &
      !$acc   grid_region_f_face,  &
      !$acc   grid_region_l_face,  &
      !$acc   grid_faces_c,  &
      !$acc   b,  &
      !$acc   flow_v_flux_n   &
      !$acc )
      do s = grid_region_f_face(reg), grid_region_l_face(reg)  ! all present
        c1 = grid_faces_c(1,s)  ! inside cell
        b(c1) = b(c1) - flow_v_flux_n(s)
      end do
      !$acc end parallel

    end if
  end do

# if T_FLOWS_DEBUG == 1
  call Grid % Save_Debug_Vtu("bp_1",               &
                             inside_name="vol_src", &
                             inside_cell=b)
# endif

  !------------------------------------------------------------------!
  !   Find the cell with the maximum volume imbalance and print it   !
  !------------------------------------------------------------------!
  max_abs_val = 0.0
  !$acc parallel loop reduction(max: max_abs_val)  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   b   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present
    max_abs_val = max(max_abs_val, abs(b(c)))
  end do
  !$acc end parallel

  ! Find maximum volume balance error over all processors
  call Global % Max_Real(max_abs_val)

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  cfl_max = 0.0
  pe_max  = 0.0

  !$acc parallel loop reduction(max: cfl_max,pe_max)  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   grid_faces_c,  &
  !$acc   visc,  &
  !$acc   dens,  &
  !$acc   flow_v_flux_n,  &
  !$acc   fc,  &
  !$acc   grid_d   &
  !$acc )
  do s = grid_region_f_face(grid_n_regions), grid_region_l_face(grid_n_regions)  ! all present
    c1 = grid_faces_c(1, s)
    c2 = grid_faces_c(2, s)

    nu_f = Face_Value(s, (visc(c1)/dens(c1)), (visc(c2)/dens(c2)))

    cfl_t   = abs(flow_v_flux_n(s)) * Flow % dt / (fc(s) * grid_d(s)**2)
    pe_t    = abs(flow_v_flux_n(s)) / fc(s) / nu_f
    cfl_max = max( cfl_max, cfl_t )
    pe_max  = max( pe_max,  pe_t  )

  end do
  !$acc end parallel

  Flow % cfl_max = cfl_max
  Flow % pe_max  = pe_max

  ! Find maximum CFL and Peclet numbers over all processors
  call Global % Max_Reals(Flow % cfl_max, Flow % pe_max)

  !-------------------------------!
  !@ Use this for REPORT_VOLUME_BALANCE somehow?
  !@ O_Print '(a,es12.3)', ' # Max. volume balance error '//  &
  !@                       'after correction: ', max_abs_val
  call Info % Iter_Fill_At(1, 5, 'dum', max_abs_val)

  call Profiler % Stop('Correct_Velocity')

  end subroutine
