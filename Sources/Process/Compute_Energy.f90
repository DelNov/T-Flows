!==============================================================================!
  subroutine Compute_Energy(flow, turb, mult, sol, dt, ini)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Field_Mod,    only: Field_Type, conductivity, capacity, density
  use Turb_Mod
  use Var_Mod,      only: Var_Type
  use Face_Mod,     only: Face_Type
  use Grid_Mod,     only: Grid_Type
  use Grad_Mod
  use Info_Mod
  use Numerics_Mod
  use Solver_Mod,   only: Solver_Type, Solver_Mod_Alias_System, Bicg, Cg, Cgs
  use Matrix_Mod,   only: Matrix_Type
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini
  real                          :: dt
!-----------------------------------[Locals]-----------------------------------! 
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: ut, vt, wt
  type(Face_Type),   pointer :: m_flux
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: n, c, s, c1, c2, exec_iter
  real                       :: a0, a12, a21
  real                       :: con_eff1, f_ex1, f_im1, tx_f1, ty_f1, tz_f1
  real                       :: con_eff2, f_ex2, f_im2, tx_f2, ty_f2, tz_f2
  real                       :: ts, pr_t1, pr_t2, pr_tf
  real                       :: ut_s, vt_s, wt_s, t_stress, con_t
  real, allocatable          :: capacity_times_density(:) 
!------------------------------------------------------------------------------!
!
!  The form of equations which are solved:
!
!     /                /                 /
!    |        dT      |                 |
!    | rho Cp -- dV   | rho u Cp T dS = | lambda  DIV T dS
!    |        dt      |                 |
!   /                /                 /
!
!
!  Dimension of the system under consideration
!
!     [A]{T} = {b}   [J/s = W]  
!
!  Dimensions of certain variables:
!
!     Cp     [J/kg K] (heat capacity)
!     lambda [W/m K] (heat conductivity)
!
!     A      [kg/s]
!     T      [K]
!     b      [kg K/s] 
!     Flux   [kg/s]
!     CT*,   [kg K/s] 
!     DT*,   [kg K/s] 
!     XT*,   [kg K/s]
! 
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  m_flux => flow % m_flux
  call Field_Mod_Alias_Momentum  (flow, u, v, w)
  call Field_Mod_Alias_Energy    (flow, t)
  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
  call Solver_Mod_Alias_System   (sol, a, b)

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0

  ! User function
  call User_Mod_Beginning_Of_Compute_Energy(flow, dt, ini)

  ! Old values (o and oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      t % oo(c) = t % o(c)
      t % o (c) = t % n(c)
    end do
  end if

  ! Gradients
  call Grad_Mod_Variable(t)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(t, capacity, m_flux % n, sol,  &
                                   t % x,                   &
                                   t % y,                   &
                                   t % z,                   &
                                   grid % dx,               &
                                   grid % dy,               &
                                   grid % dz)

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  if(turbulence_model .ne. NONE .and.  &
     turbulence_model .ne. DNS) then
  end if

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(turbulence_model .ne. LES_SMAGORINSKY    .and.  &
       turbulence_model .ne. LES_DYNAMIC        .and.  &
       turbulence_model .ne. HYBRID_LES_PRANDTL .and.  &
       turbulence_model .ne. LES_WALE           .and.  &
       turbulence_model .ne. NONE               .and.  &
       turbulence_model .ne. DNS) then
      pr_t1 = Turb_Mod_Prandtl_Number(turb, c1)
      pr_t2 = Turb_Mod_Prandtl_Number(turb, c2)
      pr_tf = grid % fw(s) * pr_t1 + (1.0-grid % fw(s)) * pr_t2
    else
      pr_tf = pr_t      
    end if

    ! Gradients on the cell face (fw corrects situation close to the wall)
    tx_f1 = grid % fw(s) * t % x(c1) + (1.0-grid % fw(s)) * t % x(c2) 
    ty_f1 = grid % fw(s) * t % y(c1) + (1.0-grid % fw(s)) * t % y(c2)
    tz_f1 = grid % fw(s) * t % z(c1) + (1.0-grid % fw(s)) * t % z(c2)
    tx_f2 = tx_f1
    ty_f2 = ty_f1
    tz_f2 = tz_f1
    if(turbulence_model .ne. NONE .and.  &
       turbulence_model .ne. DNS) then
      con_eff1 = grid % fw(s) *(conductivity+capacity*turb % vis_t(c1)/pr_tf) &
          + (1.0-grid % fw(s))*(conductivity+capacity*turb % vis_t(c2)/pr_tf)
      con_t    = grid % fw(s) *capacity*turb % vis_t(c1)/pr_tf  &
          + (1.0-grid % fw(s))*capacity*turb % vis_t(c2)/pr_tf
    else
      con_eff1 = conductivity
    end if

    con_eff2 = con_eff1

    if(turbulence_model .eq. K_EPS        .or.  &
       turbulence_model .eq. K_EPS_ZETA_F .or.  &
       turbulence_model .eq. HYBRID_LES_RANS) then
      if(c2 < 0) then
        if(Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALL .or.  &
           Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALLFL) then
          con_eff1 = turb % con_w(c1)
          con_eff2 = con_eff1
        end if
      end if
    end if

    ! Total (exact) diffusion flux
    f_ex1 = con_eff1 * (  tx_f1 * grid % sx(s)   &
                        + ty_f1 * grid % sy(s)   &
                        + tz_f1 * grid % sz(s))
    f_ex2 = con_eff2 * (  tx_f2 * grid % sx(s)   &
                        + ty_f2 * grid % sy(s)   &
                        + tz_f2 * grid % sz(s))

    ! Implicit diffusion flux
    f_im1 = con_eff1 * a % fc(s)         &
          * (  tx_f1 * grid % dx(s)      &
             + ty_f1 * grid % dy(s)      &
             + tz_f1 * grid % dz(s) )
    f_im2 = con_eff2 * a % fc(s)         &
          * (  tx_f2 * grid % dx(s)      &
             + ty_f2 * grid % dy(s)      &
             + tz_f2 * grid % dz(s) )

    ! Cross diffusion part
    t % c(c1) = t % c(c1) + f_ex1 - f_im1
    if(c2 > 0) then
      t % c(c2) = t % c(c2) - f_ex2 + f_im2
    end if

    !---------------------------!
    !                           !
    !   Turbulent heat fluxes   !
    !                           !
    !---------------------------!
    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
       turbulence_model .eq. RSM_MANCEAU_HANJALIC) then

      ! Turbulent heat fluxes according to GGDH scheme
      ! (first line is GGDH, second line is SGDH substratced 
      ut_s =  (    grid % fw(s)  * ut % n(c1)   &
           +  (1.0-grid % fw(s)) * ut % n(c2))
      vt_s =  (    grid % fw(s)  * vt % n(c1)   &
           +  (1.0-grid % fw(s)) * vt % n(c2))
      wt_s =  (    grid % fw(s)  * wt % n(c1)   &
           +  (1.0-grid % fw(s)) * wt % n(c2))
      t_stress = - (  ut_s * grid % sx(s)                  &
                    + vt_s * grid % sy(s)                  &
                    + wt_s * grid % sz(s) )                &
                    - (con_t * (  tx_f1 * grid % sx(s)     &
                                + ty_f1 * grid % sy(s)     &
                                + tz_f1 * grid % sz(s)) )

      ! Put the influence of turbulent heat fluxes explicitly in the system
      b(c1) = b(c1) + t_stress
      if(c2 > 0) then
        b(c2) = b(c2) - t_stress
      end if
    end if  ! if models are of RSM type

    ! Calculate the coefficients for the sysytem matrix
    a12 = con_eff1 * a % fc(s)
    a21 = con_eff2 * a % fc(s)

    a12 = a12  - min(m_flux % n(s), 0.0) * capacity
    a21 = a21  + max(m_flux % n(s), 0.0) * capacity

    ! Fill the system matrix
    if(c2 > 0) then
      a % val(a % dia(c1))  = a % val(a % dia(c1)) + a12
      a % val(a % dia(c2))  = a % val(a % dia(c2)) + a21
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
    else if(c2.lt.0) then
      ! Outflow is included because of the flux 
      ! corrections which also affects velocities
      if( (Var_Mod_Bnd_Cell_Type(t, c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cell_Type(t, c2) .eq. CONVECT) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * t % n(c2)
      ! In case of wallflux 
      else if(Var_Mod_Bnd_Cell_Type(t, c2) .eq. WALLFL) then
        b(c1) = b(c1) + grid % s(s) * t % q(c2)
      end if
    end if

  end do  ! through sides

  ! Cross diffusion terms are treated explicity
  do c = 1, grid % n_cells
    if(t % c(c) >= 0) then
      b(c)  = b(c) + t % c(c)
    else
      a % val(a % dia(c)) = a % val(a % dia(c))  &
                          - t % c(c) / (t % n(c) + MICRO)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  allocate(capacity_times_density(size(density)))
  capacity_times_density(:) = capacity * density(:)

  call Numerics_Mod_Inertial_Term(t, capacity_times_density, sol, dt)

  deallocate(capacity_times_density)

  !--------------------!
  !                    !
  !   User source(s)   !
  !                    !
  !--------------------!
  call User_Mod_Source(flow, t, a, b)

  !-------------------------------!
  !                               !
  !   Solve the equations for t   !
  !                               !
  !-------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(t, sol)

  ! Call linear solver to solve the equations
  call Bicg(sol,          &
            t % n,        &
            b,            &
            t % precond,  &
            t % niter,    &
            exec_iter,    &
            t % tol,      &
            t % res)

  ! Print some info on the screen
  call Info_Mod_Iter_Fill_At(1, 6, t % name, exec_iter, t % res)

  call Comm_Mod_Exchange_Real(grid, t % n)

  ! User function
  call User_Mod_End_Of_Compute_Energy(flow, dt, ini)

  end subroutine
