!==============================================================================!
  subroutine Compute_Stresses(flow, sol, dt, ini, phi, n_time_step)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equation for Re stresses for RSM.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Var_Mod,      only: Var_Type
  use Grid_Mod,     only: Grid_Type
  use Field_Mod,    only: Field_Type, density, viscosity
  use Grad_Mod
  use Info_Mod,     only: Info_Mod_Iter_Fill_At
  use Numerics_Mod, only: CENTRAL, LINEAR, PARABOLIC
  use Solver_Mod,   only: Solver_Type, Bicg, Cg, Cgs
  use Matrix_Mod,   only: Matrix_Type
  use Control_Mod
  use Work_Mod,     only: phi_x       => r_cell_01,  &
                          phi_y       => r_cell_02,  &
                          phi_z       => r_cell_03,  &
                          phi_min     => r_cell_04,  &
                          phi_max     => r_cell_05,  &
                          u1uj_phij   => r_cell_06,  &
                          u2uj_phij   => r_cell_07,  &
                          u3uj_phij   => r_cell_08,  &
                          u1uj_phij_x => r_cell_09,  &
                          u2uj_phij_y => r_cell_10,  &
                          u3uj_phij_z => r_cell_11
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
  real                      :: dt
  integer                   :: ini
  type(Var_Type)            :: phi
  integer                   :: n_time_step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w
  real,              pointer :: flux(:)
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2
  real                       :: f_ex, f_im
  real                       :: phis
  real                       :: a0, a12, a21
  real                       :: ini_res
  real                       :: vis_eff
  real                       :: phix_f, phiy_f, phiz_f
  real                       :: vis_t_f
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

  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  u    => flow % u
  v    => flow % v
  w    => flow % w
  a    => sol  % a
  b    => sol  % b % val

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
  call Grad_Mod_Component(grid, phi % n, 1, phi_x, .true.)
  call Grad_Mod_Component(grid, phi % n, 2, phi_y, .true.)
  call Grad_Mod_Component(grid, phi % n, 3, phi_z, .true.)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Turbulence(phi % adv_scheme)
  call Control_Mod_Blending_Coefficient_For_Turbulence(phi % blend)

  ! Compute phimax and phimin
  if(phi % adv_scheme .ne. CENTRAL) then
    call Calculate_Minimum_Maximum(grid, phi % n, phi_min, phi_max)
  end if

  ! New values
  do c = 1, grid % n_cells
    phi % a(c) = 0.
    phi % c(c) = 0.
  end do

  !----------------------------!
  !   Spatial Discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s) 

    ! Velocities on "orthogonal" cell centers 
    if(c2 > 0) then
      phis =        grid % f(s)  * phi % n(c1)   &
           + (1.0 - grid % f(s)) * phi % n(c2)

      ! Compute phis with desired advection scheme
      if(phi % adv_scheme .ne. CENTRAL) then
        call Advection_Scheme(flow, phis, s, phi % n, phi_min, phi_max,  &
                              phi_x, phi_y, phi_z,                       &
                              grid % dx, grid % dy, grid % dz,           &
                              phi % adv_scheme, phi % blend) 
      end if 

      ! Compute advection term
      if(c2  > 0) then
        phi % a(c1) = phi % a(c1) - flux(s)*phis
        phi % a(c2) = phi % a(c2) + flux(s)*phis
      else
        phi % a(c1) = phi % a(c1) - flux(s)*phis
      end if 

      ! Store upwinded part of the advection term in "c"
      if(flux(s)  < 0) then   ! from c2 to c1
        phi % c(c1)=phi % c(c1) - flux(s) * phi % n(c2)
        if(c2  > 0) then
          phi % c(c2)=phi % c(c2) + flux(s) * phi % n(c2)
        end if
      else
        phi % c(c1)=phi % c(c1) - flux(s) * phi % n(c1)
        if(c2  > 0) then
          phi % c(c2)=phi % c(c2) + flux(s) * phi % n(c1)
        end if
      end if
    end if     ! c2 > 0 
  end do    ! through sides

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + (phi % a(c) - phi % c(c))
  end do

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  ! Set c values back to zero
  do c = 1, grid % n_cells
    phi % c(c) = 0.
  end do

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! vis_tur is used to make diaginal element more dominant.
    ! This contribution is later substracted.
    vis_t_f = grid % fw(s) * vis_t(c1) + (1.0-grid % fw(s)) * vis_t(c2)

    vis_eff = viscosity + vis_t_f

    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      if(turbulence_model_variant .ne. STABILIZED) then
        vis_eff = 1.5*viscosity + vis_t_f
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
        u1uj_phij(c) = density * c_mu_d / phi % sigma * kin % n(c)     &
                     / max(eps % n(c), TINY)                           &
                     * (  uu % n(c) * phi_x(c)                         &
                        + uv % n(c) * phi_y(c)                         &
                        + uw % n(c) * phi_z(c))                        &
                     - viscosity * phi_x(c)

        u2uj_phij(c) = density * c_mu_d / phi % sigma * kin % n(c)     &
                     / max(eps % n(c), TINY)                           &
                     * (  uv % n(c) * phi_x(c)                         &
                        + vv % n(c) * phi_y(c)                         &
                        + vw % n(c) * phi_z(c))                        &
                     - viscosity * phi_y(c)

        u3uj_phij(c) = density * c_mu_d / phi % sigma * kin % n(c)     &
                     / max(eps % n(c), TINY)                           &
                     * (  uw % n(c) * phi_x(c)                         &
                        + vw % n(c) * phi_y(c)                         &
                        + ww % n(c) * phi_z(c))                        &
                     - viscosity * phi_z(c)
      end do
    else if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      do c = 1, grid % n_cells
        u1uj_phij(c) = density * c_mu_d / phi % sigma * t_scale(c)  &
                     * (  uu % n(c) * phi_x(c)                      &
                        + uv % n(c) * phi_y(c)                      &
                        + uw % n(c) * phi_z(c))

        u2uj_phij(c) = density * c_mu_d / phi % sigma * t_scale(c)  &
                     * (  uv % n(c) * phi_x(c)                      &
                        + vv % n(c) * phi_y(c)                      &
                        + vw % n(c) * phi_z(c))

        u3uj_phij(c) = density * c_mu_d / phi % sigma * t_scale(c)  &
                     * (  uw % n(c) * phi_x(c)                      &
                        + vw % n(c) * phi_y(c)                      &
                        + ww % n(c) * phi_z(c))
      end do
    end if

    call Grad_Mod_Component(grid, u1uj_phij, 1, u1uj_phij_x, .true.)
    call Grad_Mod_Component(grid, u2uj_phij, 2, u2uj_phij_y, .true.)
    call Grad_Mod_Component(grid, u3uj_phij, 3, u3uj_phij_z, .true.)

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

      vis_eff = (grid % fw(s) * vis_t(c1) + (1.0-grid % fw(s)) * vis_t(c2))

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

  call Control_Mod_Time_Integration_Scheme(phi % td_scheme)

  ! Two time levels; Linear interpolation
  if(phi % td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = density*grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c) = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(phi % td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = density*grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0*a0 * phi % o(c) - 0.5*a0 * phi % oo(c)
    end do
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then 
    call Grad_Mod_Component(grid, f22 % n, 1, f22 % x, .true.) ! df22/dx
    call Grad_Mod_Component(grid, f22 % n, 2, f22 % y, .true.) ! df22/dy
    call Grad_Mod_Component(grid, f22 % n, 3, f22 % z, .true.) ! df22/dz

    call Sources_Rsm_Manceau_Hanjalic(flow, sol, phi % name)
  else if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Sources_Rsm_Hanjalic_Jakirlic(flow, sol, phi % name, n_time_step)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Set under-relaxation factor then overwrite with control file if specified
  phi % urf = 1.0
  call Control_Mod_Simple_Underrelaxation_For_Turbulence(phi % urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - phi % urf)*phi % n(c)  &
         / phi % urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / phi % urf
  end do

  call Control_Mod_Tolerance_For_Turbulence_Solver(phi % tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(phi % precond)

  ! Set the number of iterations then overwrite with control file if specified
  phi % niter = 6
  call Control_Mod_Max_Iterations_For_Turbulence_Solver(phi % niter)

  call Bicg(sol,            &
            phi % n,        &
            b,              &
            phi % precond,  &
            phi % niter,    &
            phi % tol,      &
            ini_res,        &
            phi % res)

  if( phi % name .eq. 'UU' )   &
    call Info_Mod_Iter_Fill_At(3, 1, phi % name, phi % niter, phi % res)
  if( phi % name .eq. 'VV' )   &
    call Info_Mod_Iter_Fill_At(3, 2, phi % name, phi % niter, phi % res)
  if( phi % name .eq. 'WW' )   &
    call Info_Mod_Iter_Fill_At(3, 3, phi % name, phi % niter, phi % res)
  if( phi % name .eq. 'UV' )   &
    call Info_Mod_Iter_Fill_At(3, 4, phi % name, phi % niter, phi % res)
  if( phi % name .eq. 'UW' )   &
    call Info_Mod_Iter_Fill_At(3, 5, phi % name, phi % niter, phi % res)
  if( phi % name .eq. 'VW' )   &
    call Info_Mod_Iter_Fill_At(3, 6, phi % name, phi % niter, phi % res)
  if( phi % name .eq. 'EPS' )  &
    call Info_Mod_Iter_Fill_At(4, 1, phi % name, phi % niter, phi % res)

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

  end subroutine
