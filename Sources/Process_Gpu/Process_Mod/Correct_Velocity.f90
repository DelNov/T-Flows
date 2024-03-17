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
  real, contiguous, pointer :: b(:), v_flux(:)
  real, contiguous, pointer :: v_m(:), fc(:)
  type(Grid_Type),  pointer :: Grid
  real                      :: a12, b_tmp, max_abs_val
  integer                   :: c, s, c1, c2, i_cel
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Correct_Velocity')

  ! Take some aliases
  b      => Flow % Nat % b
  v_m    => Flow % Nat % M % v_m
  fc     => Flow % Nat % C % fc
  v_flux => Flow % v_flux
  Grid   => Flow % pnt_grid

  !----------------------!
  !   Correct velocity   !
  !----------------------!

  ! Units: kg m / s^2 * s / kg = m / s
  !$acc parallel loop independent
  do c = 1, grid_n_cells - grid_n_buff_cells
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
  !   Correct volume fluxes inside the domain   !
  !---------------------------------------------!

  ! Units: m * m^3 * s / kg * kg / (m s^2) = m^3 / s
  !$acc parallel loop independent
  do s = grid_reg_f_face(grid_n_regions), grid_reg_l_face(grid_n_regions)
    c1 = grid_faces_c(1, s)
    c2 = grid_faces_c(2, s)

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
  do c1 = 1, grid_n_cells - grid_n_buff_cells

    b_tmp = b(c1)
    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        b_tmp = b_tmp - v_flux(s) * merge(1,-1, c1.lt.c2)
      end if
    end do
    !$acc end loop

    ! Finish, and nullify if it is not in fluid
    b(c1) = b_tmp

  end do
  !$acc end parallel

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
  do c = 1, grid_n_cells - grid_n_buff_cells
    max_abs_val = max(max_abs_val, abs(b(c)))
  end do

  ! Find maximum volume balance error over all processors
  call Global % Max_Real(max_abs_val)

  O_Print '(a,es12.3)', ' # Max. volume balance error '//  &
                        'after correction: ', max_abs_val

  !-----------------------------------!
  !     Update the pressure field     !
  !   (hard-coded under-relaxation)   !
  !-----------------------------------!

  !$acc parallel loop independent
  do c = 1, grid_n_cells - grid_n_buff_cells
    p_n(c) = p_n(c) + 0.2 * pp_n(c)
  end do
  !$acc end parallel

  ! Update buffers for presssure over all processors
  call Grid % Exchange_Cells_Real(p_n)

  call Profiler % Stop('Correct_Velocity')

  end subroutine
