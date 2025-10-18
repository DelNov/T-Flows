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
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), fc(:), pp_x(:), pp_y(:), pp_z(:)
  real, contiguous, pointer :: visc(:), dens(:)
  real                      :: a12, b_tmp, vol_res, w1, w2
  real                      :: cfl_max, pe_max, cfl_t, pe_t, nu_f, dt
  integer                   :: c, s, c1, c2, i_cel, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Correct_Velocity')

  ! Take some aliases
  ! GPU version doesn't work if you use directly Flow % whatever_variable
  ! These aliases are really needed, not just some gimmick to shorten the code
  b    => Flow % Nat % b
  fc   => Flow % Nat % A % fc
  dens => Flow % density
  visc => Flow % viscosity
  dt   =  Flow % dt

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
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    Flow % u % n(c) = Flow % u % n(c) - pp_x(c) * Flow % v_m(c)
    Flow % v % n(c) = Flow % v % n(c) - pp_y(c) * Flow % v_m(c)
    Flow % w % n(c) = Flow % w % n(c) - pp_z(c) * Flow % v_m(c)
  end do
  !$tf-acc loop end

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
  !$tf-acc loop begin
  do s = Faces_In_Domain_And_At_Buffers()  ! all present
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    w1 = Grid % f(s)
    w2 = 1.0 - w1

    a12 = -fc(s) * (w1 * Flow % v_m(c1) + w2 * Flow % v_m(c2))

    Flow % v_flux % n(s) = Flow % v_flux % n(s)  &
                         + (Flow % pp % n(c2) - Flow % pp % n(c1)) * a12
  end do
  !$tf-acc loop end

  !-------------------------------!
  !                               !
  !   Re-compute volume sources   !
  !                               !
  !- - - - - - - - - - - - - - - -!
  !   This step is for checking   !
  !-------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()  ! all present
    b(c) = 0.0
  end do
  !$tf-acc loop end

  !---------------------------------!
  !   First consider inside faces   !
  !---------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    b_tmp = b(c1)
    do i_cel = Grid % cells_i_cells(c1),  &
               Grid % cells_n_cells(c1)

      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)

      b_tmp = b_tmp - Flow % v_flux % n(s) * merge(1,-1, c1.lt.c2)

    end do

    b(c1) = b_tmp

  end do
  !$tf-acc loop end

  !-----------------------------!
  !   Then the boundary faces   !
  !-----------------------------!

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW   .or.  &
       Grid % region % type(reg) .eq. OUTFLOW  .or.  &
       Grid % region % type(reg) .eq. PRESSURE .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)  ! inside cell
        b(c1) = b(c1) - Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if
  end do

  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    b(c) = b(c) / (Grid % vol(c) / dt)
  end do
  !$tf-acc loop end

# if T_FLOWS_DEBUG == 1
  call Grid % Save_Debug_Vtu("bp_1",               &
                             inside_name="vol_src", &
                             inside_cell=b)
# endif

  !------------------------------------------------------------------!
  !   Find the cell with the maximum volume imbalance and print it   !
  !------------------------------------------------------------------!
  vol_res = 0.0
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    vol_res = max(vol_res, abs(b(c)))
  end do
  !$tf-acc loop end
  Flow % vol_res = vol_res

  ! Find maximum volume balance error over all processors
  call Global % Max_Real(Flow % vol_res)

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  cfl_max = 0.0
  pe_max  = 0.0

  !$tf-acc loop begin
  do s = Faces_In_Domain_And_At_Buffers()  ! all present
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    w1 = Grid % f(s)
    w2 = 1.0 - w1

    nu_f = w1 * visc(c1)/dens(c1) + w2 * visc(c2)/dens(c2)

    cfl_t   = abs(Flow % v_flux % n(s)) * Flow % dt / (fc(s) * Grid % d(s)**2)
    pe_t    = abs(Flow % v_flux % n(s)) / fc(s) / nu_f
    cfl_max = max( cfl_max, cfl_t )
    pe_max  = max( pe_max,  pe_t  )

  end do
  !$tf-acc loop end

  Flow % cfl_max = cfl_max
  Flow % pe_max  = pe_max

  ! Find maximum CFL and Peclet numbers over all processors
  call Global % Max_Reals(Flow % cfl_max, Flow % pe_max)

  !-------------------------------!
  !@ Use this for REPORT_VOLUME_BALANCE somehow?
  !@ O_Print '(a,es12.3)', ' # Max. volume balance error '//  &
  !@                       'after correction: ', Flow % vol_res
  call Info % Iter_Fill_At(1, 5, 'dum', Flow % vol_res)

  call Profiler % Stop('Correct_Velocity')

  end subroutine
