!==============================================================================!
  subroutine Compute_Stress(Turb, Sol, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equation for Re stresses for RSM.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
  type(Var_Type)            :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, zeta, f22, ut, vt, wt
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  real, contiguous,  pointer :: flux(:)
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2, nc, nb
  real                       :: f_ex, f_im
  real                       :: a0, a12, a21
  real                       :: vis_eff
  real                       :: phix_f, phiy_f, phiz_f
  real                       :: vis_t_f
  real                       :: dt
  real                       :: visc_f
  real, contiguous,  pointer :: phi_x(:), phi_y(:), phi_z(:), cross(:)
  real, contiguous,  pointer :: u1uj_phij(:),   u2uj_phij(:),   u3uj_phij(:)
  real, contiguous,  pointer :: u1uj_phij_x(:), u2uj_phij_y(:), u3uj_phij_z(:)
!------------------------------------------------------------------------------!
!                                                                              !
!   The form of equations which are being solved:                              !
!                                                                              !
!     /               /                /                     /                 !
!    |     dphi      |                | mu_eff              |                  !
!    | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV             !
!    |      dt       |                |  sigma              |                  !
!   /               /                /                     /                   !
!                                                                              !
!==============================================================================!

  call Profiler % Start('Compute_Stress (without solvers)')

  call Work % Connect_Real_Cell(phi_x, phi_y, phi_z, cross,       &
                                u1uj_phij, u2uj_phij, u3uj_phij,  &
                                u1uj_phij_x, u2uj_phij_y, u3uj_phij_z)

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  nc   =  Grid % n_cells
  nb   =  Grid % n_bnd_cells
  dt   =  Flow % dt
  flux => Flow % v_flux % n
  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_Heat_Fluxes (ut, vt, wt)
  call Sol % Alias_Native       (A, b)

  ! Initialize cross diffusion sources, matrix and right hand side
  cross  (:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o) and older than old (oo)
  if(Iter % Current() .eq. 1) then
    do c = Cells_In_Domain_And_Buffers()
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
    end do
  end if

  ! Gradients
  call Flow % Grad(Grid, phi % n, phi_x(-nb:nc),  &
                                  phi_y(-nb:nc),  &
                                  phi_z(-nb:nc))

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(phi, Flow % density, flux, b)

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! vis_tur is used to make diaginal element more dominant.
    ! This contribution is later substracted.
    vis_t_f =      Grid % fw(s)  * Turb % vis_t(c1)  &
            + (1.0-Grid % fw(s)) * Turb % vis_t(c2)

    visc_f =        Grid % fw(s)  * Flow % viscosity(c1)  &
           + (1.0 - Grid % fw(s)) * Flow % viscosity(c2)

    vis_eff = visc_f + vis_t_f

    if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      if(Turb % model_variant .ne. STABILIZED) then
        vis_eff = 1.5 * visc_f + vis_t_f
      end if
    end if

    phix_f = Grid % fw(s) * phi_x(c1) + (1.0-Grid % fw(s)) * phi_x(c2)
    phiy_f = Grid % fw(s) * phi_y(c1) + (1.0-Grid % fw(s)) * phi_y(c2)
    phiz_f = Grid % fw(s) * phi_z(c1) + (1.0-Grid % fw(s)) * phi_z(c2)


    ! Total (exact) diffusive flux plus Turb. diffusion
    f_ex = vis_eff * (  phix_f * Grid % sx(s)  &
                      + phiy_f * Grid % sy(s)  &
                      + phiz_f * Grid % sz(s) ) 

    a0 = vis_eff * A % fc(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: f_coef is
    !  not corrected at interface between materials)
    f_im=( phix_f*Grid % dx(s)       &
          +phiy_f*Grid % dy(s)       &
          +phiz_f*Grid % dz(s))*a0

    ! Cross diffusion part
    cross(c1) = cross(c1) + f_ex - f_im
    if(c2  > 0) then
      cross(c2) = cross(c2) - f_ex + f_im
    end if

    ! Compute coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    ! Blend system matrix if desired to do so
    if(phi % blend_matrix) then
      a12 = a12  - min(flux(s), 0.0) * Flow % density(c1)
      a21 = a21  + max(flux(s), 0.0) * Flow % density(c2)
    end if

    ! Fill the system matrix
    if(c2  > 0) then
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % dia(c1))  = A % val(A % dia(c1))  + a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
      A % val(A % dia(c2))  = A % val(A % dia(c2))  + a21
    else if(c2  < 0) then

      ! Outflow is not included because it was causing problems     
      ! Convect is commented because for turbulent scalars convect 
      ! outflow is treated as classic outflow.
      if((Grid % Bnd_Cond_Type(c2) .eq. INFLOW).or.     &
         (Grid % Bnd_Cond_Type(c2) .eq. WALL).or.       &
!!!      (Grid % Bnd_Cond_Type(c2) .eq. CONVECT).or.    &
         (Grid % Bnd_Cond_Type(c2) .eq. WALLFL) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1) = b(c1) + a12 * phi % n(c2)
      end if
    end if

  end do  ! through faces

  !------------------------------!
  !   Turbulent diffusion term   !
  !------------------------------!
  if(phi % name .eq. 'EPS') then
    Turb % c_mu_d = 0.18
  else
    Turb % c_mu_d = 0.22
  end if

  if(Turb % model_variant .ne. STABILIZED) then
    if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      do c = Cells_In_Domain_And_Buffers()
        u1uj_phij(c) = Flow % density(c) * Turb % c_mu_d / phi % sigma &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uu % n(c) * phi_x(c)                         &
                        + uv % n(c) * phi_y(c)                         &
                        + uw % n(c) * phi_z(c))                        &
                     - Flow % viscosity(c) * phi_x(c)

        u2uj_phij(c) = Flow % density(c) * Turb % c_mu_d / phi % sigma &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uv % n(c) * phi_x(c)                         &
                        + vv % n(c) * phi_y(c)                         &
                        + vw % n(c) * phi_z(c))                        &
                     - Flow % viscosity(c) * phi_y(c)

        u3uj_phij(c) = Flow % density(c) * Turb % c_mu_d / phi % sigma &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uw % n(c) * phi_x(c)                         &
                        + vw % n(c) * phi_y(c)                         &
                        + ww % n(c) * phi_z(c))                        &
                     - Flow % viscosity(c) * phi_z(c)
      end do
    else if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
      do c = Cells_In_Domain_And_Buffers()
        u1uj_phij(c) = Flow % density(c) * Turb % c_mu_d / phi % sigma     &
                     * Turb % t_scale(c)                                   &
                     * (  uu % n(c) * phi_x(c)                             &
                        + uv % n(c) * phi_y(c)                             &
                        + uw % n(c) * phi_z(c))

        u2uj_phij(c) = Flow % density(c) * Turb % c_mu_d / phi % sigma     &
                     * Turb % t_scale(c)                                   &
                     * (  uv % n(c) * phi_x(c)                             &
                        + vv % n(c) * phi_y(c)                             &
                        + vw % n(c) * phi_z(c))

        u3uj_phij(c) = Flow % density(c) * Turb % c_mu_d / phi % sigma     &
                     * Turb % t_scale(c)                                   &
                     * (  uw % n(c) * phi_x(c)                             &
                        + vw % n(c) * phi_y(c)                             &
                        + ww % n(c) * phi_z(c))
      end do
    end if

    call Flow % Grad_Component(Grid, u1uj_phij(-nb:nc), 1, u1uj_phij_x(-nb:nc))
    call Flow % Grad_Component(Grid, u2uj_phij(-nb:nc), 2, u2uj_phij_y(-nb:nc))
    call Flow % Grad_Component(Grid, u3uj_phij(-nb:nc), 3, u3uj_phij_z(-nb:nc))

    do c = Cells_In_Domain_And_Buffers()
      b(c) = b(c) + (  u1uj_phij_x(c)  &
                     + u2uj_phij_y(c)  &
                     + u3uj_phij_z(c) ) * Grid % vol(c)
    end do
  end if

  !------------------------------------------------------------------!
  !   Here we clean up transport equation from the false diffusion   !
  !------------------------------------------------------------------!
  if(Turb % model_variant .ne. STABILIZED) then
    do s = 1, Grid % n_faces

      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      vis_eff = (Grid % fw(s)      * Turb % vis_t(c1)  &
              + (1.0-Grid % fw(s)) * Turb % vis_t(c2))

      phix_f = Grid % fw(s) * phi_x(c1) + (1.0-Grid % fw(s)) * phi_x(c2)
      phiy_f = Grid % fw(s) * phi_y(c1) + (1.0-Grid % fw(s)) * phi_y(c2)
      phiz_f = Grid % fw(s) * phi_z(c1) + (1.0-Grid % fw(s)) * phi_z(c2)
      f_ex = vis_eff * (  phix_f * Grid % sx(s)  &
                        + phiy_f * Grid % sy(s)  &
                        + phiz_f * Grid % sz(s))
      a0 = vis_eff * A % fc(s)
      f_im = (   phix_f * Grid % dx(s)        &
               + phiy_f * Grid % dy(s)        &
               + phiz_f * Grid % dz(s)) * a0

      b(c1) = b(c1)                                             &
             - vis_eff * (phi % n(c2) - phi%n(c1)) * A % fc(s)  &
             - f_ex + f_im
      if(c2  > 0) then
        b(c2) = b(c2)                                            &
              + vis_eff * (phi % n(c2) - phi%n(c1)) * A % fc(s)  &
              + f_ex - f_im
      end if
    end do
  end if

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()
    b(c) = b(c) + cross(c)
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  call Numerics_Mod_Inertial_Term(phi, Flow % density, a, b, dt)

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !                                     !
  !-------------------------------------!
  if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Flow % Grad_Variable(f22)

    call Turb % Src_Rsm_Manceau_Hanjalic(Sol, phi % name)
  else if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb % Src_Rsm_Hanjalic_Jakirlic(Sol, phi % name)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, a, b)

  call Profiler % Start(String % First_Upper(phi % solver)  //  &
                        ' (solver for turbulence)')

  ! Call linear solver to solve the equations
  call Sol % Run(A, phi, b)

  call Profiler % Stop(String % First_Upper(phi % solver)  //  &
                       ' (solver for turbulence)')

  ! Print info on the screen
  if( phi % name .eq. 'UU' )   &
    call Info % Iter_Fill_At(3, 1, phi % name, phi % res, phi % niter)
  if( phi % name .eq. 'VV' )   &
    call Info % Iter_Fill_At(3, 2, phi % name, phi % res, phi % niter)
  if( phi % name .eq. 'WW' )   &
    call Info % Iter_Fill_At(3, 3, phi % name, phi % res, phi % niter)
  if( phi % name .eq. 'UV' )   &
    call Info % Iter_Fill_At(3, 4, phi % name, phi % res, phi % niter)
  if( phi % name .eq. 'UW' )   &
    call Info % Iter_Fill_At(3, 5, phi % name, phi % res, phi % niter)
  if( phi % name .eq. 'VW' )   &
    call Info % Iter_Fill_At(3, 6, phi % name, phi % res, phi % niter)
  if( phi % name .eq. 'EPS' )  &
    call Info % Iter_Fill_At(4, 1, phi % name, phi % res, phi % niter)

  if(phi % name .eq. 'EPS') then
    do c= 1, Grid % n_cells
      phi % n(c) = phi % n(c)
     if( phi % n(c) < 0.) then
       phi % n(c) = phi % o(c)
     end if
    end do
  end if

  if(phi % name .eq. 'UU' .or.  &
     phi % name .eq. 'VV' .or.  &
     phi % name .eq. 'WW') then
    do c = Cells_In_Domain_And_Buffers()
      phi % n(c) = phi % n(c)
      if(phi % n(c) < 0.) then
        phi % n(c) = phi % o(c)
      end if
    end do
  end if

  call Flow % Grad_Variable(phi)

  call Work % Disconnect_Real_Cell(phi_x, phi_y, phi_z, cross,       &
                                   u1uj_phij, u2uj_phij, u3uj_phij,  &
                                   u1uj_phij_x, u2uj_phij_y, u3uj_phij_z)

  call Profiler % Stop('Compute_Stress (without solvers)')

  end subroutine
