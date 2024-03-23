!==============================================================================!
  subroutine Correct_Velocity(Process, Flow, Grid)
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
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), v_flux(:)
  real, contiguous, pointer :: v_m(:), fc(:)
  real                      :: a12, b_tmp, max_abs_val
  integer                   :: c, s, c1, c2, i_cel, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Correct_Velocity')

  ! Take some aliases
  b      => Flow % Nat % b
  fc     => Flow % Nat % C % fc
  v_flux => Flow % v_flux
  v_m    => Flow % v_m

  !----------------------!
  !                      !
  !   Correct velocity   !
  !                      !
  !----------------------!

  ! Units: kg m / s^2 * s / kg = m / s
  !$acc parallel loop independent
  do c = Cells_In_Domain()
    u_n(c) = u_n(c) - pp_x(c) * v_m(c)
    v_n(c) = v_n(c) - pp_y(c) * v_m(c)
    w_n(c) = w_n(c) - pp_z(c) * v_m(c)
  end do
  !$acc end parallel

  ! Update buffers for velocities over all processors
  call Grid % Exchange_Cells_Real(u_n)
  call Grid % Exchange_Cells_Real(v_n)
  call Grid % Exchange_Cells_Real(w_n)

  !---------------------------------------------!
  !                                             !
  !   Correct volume fluxes inside the domain   !
  !                                             !
  !---------------------------------------------!

  ! Units: m * m^3 * s / kg * kg / (m s^2) = m^3 / s
  !$acc parallel loop independent
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    a12 = -fc(s) * 0.5 * (v_m(c1) + v_m(c2))

    v_flux(s) = v_flux(s) + (pp_n(c2) - pp_n(c1)) * a12
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

  !$acc parallel loop independent
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

      !$acc parallel loop
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)  ! inside cell
        b(c1) = b(c1) - v_flux(s)
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
  !$acc parallel loop reduction(max:max_abs_val)
  do c = Cells_In_Domain()
    max_abs_val = max(max_abs_val, abs(b(c)))
  end do

  ! Find maximum volume balance error over all processors
  call Global % Max_Real(max_abs_val)

  !@ Use this for REPORT_VOLUME_BALANCE somehow?
  !@ O_Print '(a,es12.3)', ' # Max. volume balance error '//  &
  !@                       'after correction: ', max_abs_val
  call Info % Iter_Fill_At(1, 5, 'dum', max_abs_val)

  call Profiler % Stop('Correct_Velocity')

  end subroutine
