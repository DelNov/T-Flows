!==============================================================================!
  subroutine Correct_Velocity(Process, Flow, Vof, Sol)
!------------------------------------------------------------------------------!
!>  The subroutine Correct_Velocity in T-Flows is dedicated to adjusting the
!>  velocities and volume flow rates on cell faces. This is a crucial step
!>  in ensuring that the computed flow field conserves volume.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization and setup: Begins by setting up necessary variables and   !
!     pointers, including those for the grid, flow field, and solver matrices. !
!     It also initializes the profiler for performance tracking.               !
!   * Velocity and flux correction: Implements corrections to the velocity     !
!     fields (u, v, w) using the pressure corrections calculated in the        !
!     pressure solver. These corrections are essential to ensure a divergence- !
!     free velocity field, satisfying the continuity equation.                 !
!   * Flux correction at faces: Adjusts the volume flow rates at the           !
!     internal faces of cells, based on the pressure corrections, to maintain  !
!     the consistency and accuracy of the flow field.                          !
!   * Error calculation: Calculates the maximum volume error with the newly    !
!     corrected fluxes, providing a measure of the divergence in the flow.     !
!   * CFL and Peclet number calculation: Computes the CFL and Peclet numbers   !
!     for the flow, which are critical for assessing the stability and         !
!     convection dominance in the numerical simulation.                        !
!   * User-defined functions: Calls user-defined functions at specific points  !
!     for additional customization and integration of specific behaviors into  !
!     the velocity correction process.                                         !
!   * Performance monitoring: Continuously monitors the performance of the     !
!     subroutine, aiding in optimization and analysis of the simulation.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Vof_Type),      target :: Vof      !! VOF object
  type(Solver_Type),   target :: Sol      !! solver object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: A               ! pressure matrix
  type(Matrix_Type), pointer :: M               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: cfl_t, pe_t, dens_f, visc_f, dt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Correct_Velocity')

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  v_flux => Flow % v_flux
  p      => Flow % p
  pp     => Flow % pp
  dt     =  Flow % dt
  A      => Sol % Nat % A
  M      => Sol % Nat % M
  b      => Sol % Nat % b % val
  call Flow % Alias_Momentum(u, v, w)

  ! User function
  call User_Mod_Beginning_Of_Correct_Velocity(Flow, Vof, Sol)

  !-----------------------------------------!
  !   Correct velocities and fluxes with    !
  !    periodic part of the pressure to     !
  !    obtain divergence free velocity      !
  !-----------------------------------------!

  ! Normal correction
  do c = 1, Grid % n_cells
    u % n(c) = u % n(c) - pp % x(c) * Grid % vol(c) / M % sav(c)
    v % n(c) = v % n(c) - pp % y(c) * Grid % vol(c) / M % sav(c)
    w % n(c) = w % n(c) - pp % z(c) * Grid % vol(c) / M % sav(c)
  end do

  !----------------------------------------------------------------!
  !   Look at the following equation and you will understand why   !
  !   is the matrix for pressure corrections in SIMPLE algorithm   !
  !   formed from the coefficients of the velocity matrix.         !
  !   pp      [kg/ms^2]                                            !
  !   A % val [m^4s/kg]                                            !
  !----------------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then
                                     !<--- this is correction --->!
      v_flux % n(s) = v_flux % n(s) + ( pp % n(c2) - pp % n(c1) )   &
                                       * A % val(A % pos(1,s))
    end if
  end do

  !------------------------------------!
  !   Calculate the max volume error   !
  !   with the new corrected fluxes    !
  !------------------------------------!
  b(:) = 0.0

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    b(c1) = b(c1) - v_flux % n(s)
    if(c2 > 0) then
      b(c2) = b(c2) + v_flux % n(s)
    end if
  end do

  ! In case of mass transfer, also add
  ! volume change due to mass transfer
  call Vof % Mass_Transfer_Pressure_Source(b)

  do c = 1, Grid % n_cells
    b(c) = b(c) / (Grid % vol(c) / dt)
  end do

  Flow % vol_res = 0.0
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    Flow % vol_res = max(Flow % vol_res, abs(b(c)))
  end do
  call Global % Max_Real(Flow % vol_res)

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  Flow % cfl_max = 0.0
  Flow % pe_max  = 0.0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    dens_f =        Grid % fw(s)  * Flow % density  (c1)   &
           + (1.0 - Grid % fw(s)) * Flow % density  (c2)
    visc_f =        Grid % fw(s)  * Flow % viscosity(c1)   &
           + (1.0 - Grid % fw(s)) * Flow % viscosity(c2)
    if(c2 > 0) then
      cfl_t = abs( dt * v_flux % n(s) /          &
                   ( A % fc(s) *                 &
                   (  Grid % dx(s)*Grid % dx(s)  &
                    + Grid % dy(s)*Grid % dy(s)  &
                    + Grid % dz(s)*Grid % dz(s)) ) )
      pe_t    = abs( v_flux % n(s) / A % fc(s) / (visc_f / dens_f + TINY) )
      Flow % cfl_max = max( Flow % cfl_max, cfl_t )
      Flow % pe_max  = max( Flow % pe_max,  pe_t  )
    end if
  end do
  call Global % Max_Real(Flow % cfl_max)
  call Global % Max_Real(Flow % pe_max)

  if(Flow % p_m_coupling == SIMPLE) then
    call Info % Iter_Fill_At(1, 5, 'dum', Flow % vol_res)
  else
    if(Flow % i_corr == Flow % n_piso_corrections) then
      call Info % Iter_Fill_At(1, 5, 'dum', Flow % vol_res)
    end if
  end if

  ! User function
  call User_Mod_End_Of_Correct_Velocity(Flow, Vof, Sol)

  call Profiler % Stop('Correct_Velocity')

  end subroutine
