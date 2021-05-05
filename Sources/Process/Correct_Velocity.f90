!==============================================================================!
  subroutine Correct_Velocity(flow, turb, mult, sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Corrects the velocities, and mass (or volume) fluxes on cell faces.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: curr_dt
  integer, intent(in)           :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: a               ! pressure matrix
  type(Matrix_Type), pointer :: m               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: c, c1, c2, s, i_fac
  real                       :: cfl_t, pe_t, dens_f, visc_f, dt, wght
!==============================================================================!

  call Cpu_Timer_Mod_Start('Correct_Velocity')

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  v_flux => flow % v_flux
  p      => flow % p
  pp     => flow % pp
  dt     =  flow % dt
  a      => sol % a
  m      => sol % m
  b      => sol % b % val

  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  call User_Mod_Beginning_Of_Correct_Velocity(flow, mult, sol, curr_dt, ini)

  !-----------------------------------------!
  !   Correct velocities and fluxes with    !
  !    periodic part of the pressure to     !
  !    obtain divergence free velocity      !
  !-----------------------------------------!

  ! Normal correction
  do c = 1, grid % n_cells
    u % n(c) = u % n(c) - pp % x(c) * grid % vol(c) / M % sav(c)
    v % n(c) = v % n(c) - pp % y(c) * grid % vol(c) / M % sav(c)
    w % n(c) = w % n(c) - pp % z(c) * grid % vol(c) / M % sav(c)
  end do

  !----------------------------------------------------------------!
  !   Look at the following equation and you will understand why   !
  !   is the matrix for pressure corrections in SIMPLE algorithm   !
  !   formed from the coefficients of the velocity matrix.         !
  !   pp      [kg/ms^2]                                            !
  !   A % val [m^4s/kg]                                            !
  !----------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then
                                     !<--- this is correction --->!
      v_flux % n(s) = v_flux % n(s) + ( pp % n(c2) - pp % n(c1) )   &
                                       * A % val(A % pos(1,s))
    end if
  end do

  !-------------------------------------!
  !    Calculate the max mass error     !
  !   with the new (corrected) fluxes   !
  !-------------------------------------!
  do c = 1, grid % n_cells
    b(c) = 0.0
  end do

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    b(c1) = b(c1) - v_flux % n(s)
    if(c2 > 0) then
      b(c2) = b(c2) + v_flux % n(s)
    end if
  end do

  do c = 1, grid % n_cells
    b(c) = b(c) / (grid % vol(c) / dt)
  end do

  flow % vol_res = 0.0
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    flow % vol_res = max(flow % vol_res, abs(b(c)))
  end do
  call Comm_Mod_Global_Max_Real(flow % vol_res)

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  flow % cfl_max = 0.0
  flow % pe_max  = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    dens_f =        grid % fw(s)  * flow % density  (c1)   &
           + (1.0 - grid % fw(s)) * flow % density  (c2)
    visc_f =        grid % fw(s)  * flow % viscosity(c1)   &
           + (1.0 - grid % fw(s)) * flow % viscosity(c2)
    if(c2 > 0) then
      cfl_t = abs( dt * v_flux % n(s) /          &
                   ( A % fc(s) *                 &
                   (  grid % dx(s)*grid % dx(s)  &
                    + grid % dy(s)*grid % dy(s)  &
                    + grid % dz(s)*grid % dz(s)) ) )
      pe_t    = abs( v_flux % n(s) / A % fc(s) / (visc_f / dens_f + TINY) )
      flow % cfl_max = max( flow % cfl_max, cfl_t )
      flow % pe_max  = max( flow % pe_max,  pe_t  )
    end if
  end do
  call Comm_Mod_Global_Max_Real(flow % cfl_max)
  call Comm_Mod_Global_Max_Real(flow % pe_max)

  if (flow % p_m_coupling == SIMPLE) then
    call Info_Mod_Iter_Fill_At(1, 5, 'dum', -1, flow % vol_res)
  else
    if (flow % i_corr == flow % n_piso_corrections) then
      call Info_Mod_Iter_Fill_At(1, 5, 'dum', -1, flow % vol_res)
    end if
  end if

  ! User function
  call User_Mod_End_Of_Correct_Velocity(flow, mult, sol, curr_dt, ini)

  call Cpu_Timer_Mod_Stop('Correct_Velocity')

  end subroutine
