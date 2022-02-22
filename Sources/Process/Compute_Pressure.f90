!==============================================================================!
  subroutine Compute_Pressure(Flow, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the SIMPLE method.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: A               ! pressure matrix
  type(Matrix_Type), pointer :: M               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  integer, save              :: total_cells
  real                       :: total_source
  real                       :: p_max, p_min, p_nor, p_nor_c, dt, a12, fs
  character(SL)              :: solver
!==============================================================================!
!
!   The form of equations which I am solving:
!
!      /           /
!     |           |
!     | u dS = dt | GRAD pp dS
!     |           |
!    /           /
!
!   Dimension of the system under consideration
!
!     [App] {pp} = {bpp}     [m^3/s]
!
!   Dimensions of certain variables
!
!     app                   [m^4s/kg]
!     pp                    [kg/ms^2]
!     p % x, p % y, p % z   [kg/(m^2 s^2)]
!     px_f, py_f, pz_f      [kg/(m^2 s^2)]
!     b                     [m^3/s]
!     v_flux                [m^3/s]
!
!==============================================================================!

  call Cpu_Timer % Start('Compute_Pressure (without solvers)')

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
  call User_Mod_Beginning_Of_Compute_Pressure(Flow, Vof, Sol, curr_dt, ini)

  !--------------------------------------------------!
  !   Find the value for normalization of pressure   !
  !--------------------------------------------------!

  ! From control file
  call Control_Mod_Normalization_For_Pressure_Solver(p_nor_c)

  ! Calculate pressure magnitude for normalization of pressure solution
  p_max = -HUGE
  p_min = +HUGE
  do c = 1, Grid % n_cells
    p_max = max(p_max, p % n(c))
    p_min = min(p_min, p % n(c))
  end do
  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  ! Normalize pressure with the maximum of pressure difference, 
  ! value defined in control file and pressure drops.
  p_nor = max( (p_max-p_min), p_nor_c, abs(bulk % p_drop_x),  &
                                       abs(bulk % p_drop_y),  &
                                       abs(bulk % p_drop_z) )

  ! Initialize matrix and right hand side
  b       = 0.0
  A % val = 0.0

  !-----------------------------------------!
  !   Initialize the pressure corrections   !
  !    (Convergence is faster with this)    !
  !-----------------------------------------!
  pp % n = 0.0

  !----------------------------------!
  !   Correct fluxes at boundaries   !
  !----------------------------------!
  call Balance_Volume(Flow, Vof)

  !---------------------------------------!
  !   Compute volume fluxes at internal   !
  !    faces with Rhie and Chow method    !
  !---------------------------------------!
  call Rhie_And_Chow(Flow, Vof, Sol)

  !------------------------------------------!
  !   Update fluxes at boundaries and fill   !
  !   up source term for pressure equation   !
  !------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    fs = Grid % f(s)

    ! Internal fluxes fixed with Rhie and Chow method, just update source
    if(c2 > 0) then

      b(c1) = b(c1) - v_flux % n(s)
      b(c2) = b(c2) + v_flux % n(s)

    ! Side is on the boundary
    else

      if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW) then

        v_flux % n(s) = ( u % n(c2) * Grid % sx(s)     &
                        + v % n(c2) * Grid % sy(s)     &
                        + w % n(c2) * Grid % sz(s) )

        b(c1) = b(c1) - v_flux % n(s)

      else if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW  .or.  &
              Grid % Bnd_Cond_Type(c2) .eq. PRESSURE .or.  &
              Grid % Bnd_Cond_Type(c2) .eq. CONVECT) then

        v_flux % n(s) = ( u % n(c2) * Grid % sx(s)     &
                        + v % n(c2) * Grid % sy(s)     &
                        + w % n(c2) * Grid % sz(s) )

        b(c1) = b(c1) - v_flux % n(s)

        a12 = A % fc(s) * Grid % vol(c1) / M % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12

      else  ! it is SYMMETRY
        v_flux % n(s) = 0.0
      end if
    end if

  end do

  !-------------------------------------------------------------------------!
  !   In case of mass transfer, add addtional source to pressure equation   !
  !-------------------------------------------------------------------------!
  call Vof % Mass_Transfer_Pressure_Source(b)

  !----------------------------------------!
  !   Balance the source over processors   !
  !----------------------------------------!
  ! Since the introduction of face-based body forces, a very slight
  ! imbalance in pressure source can occur globally (over all processors)
  ! and we are talking about the values of the order 1e-18 here.  This
  ! happens because of the round-off errors associated with real arithmetics
  ! performed in CPUs, meaning that: A + B = B + A + epsilon.  If forces
  ! are computed over faces, they are slightly different in different
  ! processors due to these round-off errors, although values in cells
  ! surrounding them are the same after exchanging the buffers.  In order
  ! to fix it, the balancing procedure which follows is introduced.
  ! However, the fix should not be applied if domain has pressure outlet.
  !
  ! Update on July 17, 2021: I have some reservations about this part, since
  ! there was another bug fix when computing fluxes in the meanwhile (check:
  ! 90f77a1c8bd4ca05330a4435ed6321782ef00199).  This balancing also caused a
  ! bug when loading backup file (also check "Initialize_Variables" and 
  ! "Backup_Mod/Load and Backup_Mod/Save" procedures)
  if( .not. Flow % has_pressure_outlet) then
    if(total_cells .eq. 0) then  ! wasn't set yet
      total_cells = Grid % n_cells - Grid % comm % n_buff_cells
      call Comm_Mod_Global_Sum_Int(total_cells)
    end if
    total_source = sum(b(1:(Grid % n_cells - Grid % comm % n_buff_cells)))
    call Comm_Mod_Global_Sum_Real(total_source)
    b(:) = b(:) - total_source/real(total_cells)
  end if

  ! Get solver
  call Control_Mod_Solver_For_Pressure(solver)

  call Cpu_Timer % Start('Linear_Solver_For_Pressure')

  ! Call linear solver
  call Sol % Run(pp % solver,     &
                 pp % prec,       &
                 pp % prec_opts,  &
                 A,               &
                 pp % n,          &
                 b,               &
                 pp % mniter,     &  ! max number of iterations
                 pp % eniter,     &  ! executed number of iterations
                 pp % tol,        &  ! tolerance
                 pp % res,        &  ! final residual
                 norm = p_nor)       ! number for normalisation

  call Cpu_Timer % Stop('Linear_Solver_For_Pressure')

  if (Flow % p_m_coupling == SIMPLE) then
    call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
  else
    if (Flow % i_corr == Flow % n_piso_corrections) then
      call Info_Mod_Iter_Fill_At(1, 4, pp % name, pp % eniter, pp % res)
    end if
  end if

  call Flow % Grad_Pressure_Correction(pp)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    p % n(c) =  p % n(c) + pp % urf * pp % n(c)
  end do

  call Flow % Grad_Pressure(p)

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  p_max  = maxval(p % n(1:Grid % n_cells))
  p_min  = minval(p % n(1:Grid % n_cells))

  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  p % n(:) = p % n(:) - 0.5*(p_max+p_min)

  ! User function
  call User_Mod_End_Of_Compute_Pressure(Flow, Vof, Sol, curr_dt, ini)

  call Cpu_Timer % Stop('Compute_Pressure (without solvers)')

  end subroutine
