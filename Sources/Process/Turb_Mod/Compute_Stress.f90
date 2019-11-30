!==============================================================================!
  subroutine Turb_Mod_Compute_Stress(turb, sol, ini, phi, n_time_step)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equation for Re stresses for RSM.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_x       => r_cell_01,  &
                      phi_y       => r_cell_02,  &
                      phi_z       => r_cell_03,  &
                      u1uj_phij   => r_cell_06,  &
                      u2uj_phij   => r_cell_07,  &
                      u3uj_phij   => r_cell_08,  &
                      u1uj_phij_x => r_cell_09,  &
                      u2uj_phij_y => r_cell_10,  &
                      u3uj_phij_z => r_cell_11
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
  integer                   :: ini
  type(Var_Type)            :: phi
  integer                   :: n_time_step
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, zeta, f22, ut, vt, wt
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  real,              pointer :: flux(:)
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, exec_iter
  real                       :: f_ex, f_im
  real                       :: phis
  real                       :: a0, a12, a21
  real                       :: vis_eff
  real                       :: phix_f, phiy_f, phiz_f
  real                       :: vis_t_f
  real                       :: dt
  real                       :: visc_const, dens_const
!==============================================================================!
!                                                                              !
!   The form of equations which are being solved:                              !
!                                                                              !
!     /               /                /                     /                 !
!    |     dphi      |                | mu_eff              |                  !
!    | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV             !
!    |      dt       |                |  sigma              |                  !
!   /               /                /                     /                   !
!                                                                              !
!------------------------------------------------------------------------------!

  call Cpu_Timer_Mod_Start('Compute_Turbulence (without solvers)')

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  dt   =  flow % dt
  flux => flow % m_flux % n
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_Heat_Fluxes (turb, ut, vt, wt)
  call Solver_Mod_Alias_System    (sol, a, b)

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
    end do
  end if

  ! Gradients
  call Field_Mod_Grad_Component(flow, phi % n, 1, phi_x)
  call Field_Mod_Grad_Component(flow, phi % n, 2, phi_y)
  call Field_Mod_Grad_Component(flow, phi % n, 3, phi_z)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(phi, 1.0, flux, sol,  &
                                   phi % x,              &
                                   phi % y,              &
                                   phi % z,              &
                                   grid % dx,            &
                                   grid % dy,            &
                                   grid % dz)

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! vis_tur is used to make diaginal element more dominant.
    ! This contribution is later substracted.
    vis_t_f = grid % fw(s)       * turb % vis_t(c1)  &
            + (1.0-grid % fw(s)) * turb % vis_t(c2)

    visc_const = grid % f(s)         * flow % viscosity(c1)  &
               + (1.0 - grid % f(s)) * flow % viscosity(c2)

    vis_eff = visc_const + vis_t_f

    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      if(turbulence_model_variant .ne. STABILIZED) then
        vis_eff = 1.5*visc_const + vis_t_f
      end if
    end if

    phix_f = grid % fw(s) * phi_x(c1) + (1.0-grid % fw(s)) * phi_x(c2)
    phiy_f = grid % fw(s) * phi_y(c1) + (1.0-grid % fw(s)) * phi_y(c2)
    phiz_f = grid % fw(s) * phi_z(c1) + (1.0-grid % fw(s)) * phi_z(c2)


    ! Total (exact) diffusive flux plus turb. diffusion
    f_ex = vis_eff * (  phix_f * grid % sx(s)  &
                      + phiy_f * grid % sy(s)  &
                      + phiz_f * grid % sz(s) ) 

    a0 = vis_eff * a % fc(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: f_coef is
    !  not corrected at interface between materials)
    f_im=( phix_f*grid % dx(s)                      &
         +phiy_f*grid % dy(s)                      &
         +phiz_f*grid % dz(s))*a0

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex - f_im
    if(c2  > 0) then
      phi % c(c2) = phi % c(c2) - f_ex + f_im
    end if

    ! Compute coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    a12 = a12  - min(flux(s), 0.)
    a21 = a21  + max(flux(s), 0.)

    ! Fill the system matrix
    if(c2  > 0) then
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % dia(c1))  = a % val(a % dia(c1))  + a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
      a % val(a % dia(c2))  = a % val(a % dia(c2))  + a21
    else if(c2  < 0) then

      ! Outflow is not included because it was causing problems     
      ! Convect is commented because for turbulent scalars convect 
      ! outflow is treated as classic outflow.
      if((Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW).or.     &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL).or.       &
!!!      (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT).or.    &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1) = b(c1) + a12 * phi % n(c2)
      end if
    end if

  end do  ! through faces

  !------------------------------!
  !   Turbulent diffusion term   !
  !------------------------------!
  if(phi % name .eq. 'EPS') then
    c_mu_d = 0.18
  else
    c_mu_d = 0.22
  end if

  if(turbulence_model_variant .ne. STABILIZED) then
    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      do c = 1, grid % n_cells
        u1uj_phij(c) = flow % density(c) * c_mu_d / phi % sigma        &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uu % n(c) * phi_x(c)                         &
                        + uv % n(c) * phi_y(c)                         &
                        + uw % n(c) * phi_z(c))                        &
                     - flow % viscosity(c) * phi_x(c)

        u2uj_phij(c) = flow % density(c) * c_mu_d / phi % sigma        &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uv % n(c) * phi_x(c)                         &
                        + vv % n(c) * phi_y(c)                         &
                        + vw % n(c) * phi_z(c))                        &
                     - flow % viscosity(c) * phi_y(c)

        u3uj_phij(c) = flow % density(c) * c_mu_d / phi % sigma        &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uw % n(c) * phi_x(c)                         &
                        + vw % n(c) * phi_y(c)                         &
                        + ww % n(c) * phi_z(c))                        &
                     - flow % viscosity(c) * phi_z(c)
      end do
    else if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      do c = 1, grid % n_cells
        u1uj_phij(c) = flow % density(c) * c_mu_d / phi % sigma            &
                     * turb % t_scale(c)                                   &
                     * (  uu % n(c) * phi_x(c)                             &
                        + uv % n(c) * phi_y(c)                             &
                        + uw % n(c) * phi_z(c))

        u2uj_phij(c) = flow % density(c) * c_mu_d / phi % sigma            &
                     * turb % t_scale(c)                                   &
                     * (  uv % n(c) * phi_x(c)                             &
                        + vv % n(c) * phi_y(c)                             &
                        + vw % n(c) * phi_z(c))

        u3uj_phij(c) = flow % density(c) * c_mu_d / phi % sigma            &
                     * turb % t_scale(c)                                   &
                     * (  uw % n(c) * phi_x(c)                             &
                        + vw % n(c) * phi_y(c)                             &
                        + ww % n(c) * phi_z(c))
      end do
    end if

    call Field_Mod_Grad_Component(flow, u1uj_phij, 1, u1uj_phij_x)
    call Field_Mod_Grad_Component(flow, u2uj_phij, 2, u2uj_phij_y)
    call Field_Mod_Grad_Component(flow, u3uj_phij, 3, u3uj_phij_z)

    do c = 1, grid % n_cells
      b(c) = b(c) + (  u1uj_phij_x(c)  &
                     + u2uj_phij_y(c)  &
                     + u3uj_phij_z(c) ) * grid % vol(c)
    end do
  end if

  !------------------------------------------------------------------!
  !   Here we clean up transport equation from the false diffusion   !
  !------------------------------------------------------------------!
  if(turbulence_model_variant .ne. STABILIZED) then
    do s = 1, grid % n_faces

      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      vis_eff = (grid % fw(s)      * turb % vis_t(c1)  &
              + (1.0-grid % fw(s)) * turb % vis_t(c2))

      phix_f = grid % fw(s) * phi_x(c1) + (1.0-grid % fw(s)) * phi_x(c2)
      phiy_f = grid % fw(s) * phi_y(c1) + (1.0-grid % fw(s)) * phi_y(c2)
      phiz_f = grid % fw(s) * phi_z(c1) + (1.0-grid % fw(s)) * phi_z(c2)
      f_ex = vis_eff * (  phix_f * grid % sx(s)  &
                        + phiy_f * grid % sy(s)  &
                        + phiz_f * grid % sz(s))
      a0 = vis_eff * a % fc(s)
      f_im = (   phix_f * grid % dx(s)        &
               + phiy_f * grid % dy(s)        &
               + phiz_f * grid % dz(s)) * a0

      b(c1) = b(c1)                                             &
             - vis_eff * (phi % n(c2) - phi%n(c1)) * a % fc(s)  &
             - f_ex + f_im
      if(c2  > 0) then
        b(c2) = b(c2)                                            &
              + vis_eff * (phi % n(c2) - phi%n(c1)) * a % fc(s)  &
              + f_ex - f_im
      end if
    end do
  end if

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + phi % c(c)
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  call Numerics_Mod_Inertial_Term(phi, flow % density, sol, dt)

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !                                     !
  !-------------------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Field_Mod_Grad_Variable(flow, f22)

    call Turb_Mod_Src_Rsm_Manceau_Hanjalic(turb, sol, phi % name)
  else if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb_Mod_Src_Rsm_Hanjalic_Jakirlic(turb, sol, phi % name, n_time_step)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, sol)

  ! Call linear solver to solve the equations
  call Cpu_Timer_Mod_Start('Linear_Solver_For_Turbulence')
  call Bicg(sol,            &
            phi % n,        &
            b,              &
            phi % precond,  &
            phi % niter,    &
            exec_iter,      &
            phi % tol,      &
            phi % res)
  call Cpu_Timer_Mod_Stop('Linear_Solver_For_Turbulence')

  ! Print info on the screen
  if( phi % name .eq. 'UU' )   &
    call Info_Mod_Iter_Fill_At(3, 1, phi % name, exec_iter, phi % res)
  if( phi % name .eq. 'VV' )   &
    call Info_Mod_Iter_Fill_At(3, 2, phi % name, exec_iter, phi % res)
  if( phi % name .eq. 'WW' )   &
    call Info_Mod_Iter_Fill_At(3, 3, phi % name, exec_iter, phi % res)
  if( phi % name .eq. 'UV' )   &
    call Info_Mod_Iter_Fill_At(3, 4, phi % name, exec_iter, phi % res)
  if( phi % name .eq. 'UW' )   &
    call Info_Mod_Iter_Fill_At(3, 5, phi % name, exec_iter, phi % res)
  if( phi % name .eq. 'VW' )   &
    call Info_Mod_Iter_Fill_At(3, 6, phi % name, exec_iter, phi % res)
  if( phi % name .eq. 'EPS' )  &
    call Info_Mod_Iter_Fill_At(4, 1, phi % name, exec_iter, phi % res)

  if(phi % name .eq. 'EPS') then
    do c= 1, grid % n_cells
      phi % n(c) = phi % n(c)
     if( phi % n(c) < 0.) then
       phi % n(c) = phi % o(c)
     end if
    end do
  end if

  if(phi % name .eq. 'UU' .or.  &
     phi % name .eq. 'VV' .or.  &
     phi % name .eq. 'WW') then
    do c = 1, grid % n_cells
      phi % n(c) = phi % n(c)
      if(phi % n(c) < 0.) then
        phi % n(c) = phi % o(c)
      end if
    end do
  end if

  call Comm_Mod_Exchange_Real(grid, phi % n)

  call Cpu_Timer_Mod_Stop('Compute_Turbulence (without solvers)')

  end subroutine
