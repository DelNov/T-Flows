!==============================================================================!
  subroutine Compute_Stresses(grid, dt, ini, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equation for Re stresses for RSM.         !
!   'EBM' and 'HJ' are calling this subroutine.                                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Comm_Mod
  use Var_Mod
  use Grid_Mod
  use Grad_Mod
  use Info_Mod
  use Numerics_Mod
  use Solvers_Mod, only: Bicg, Cg, Cgs
  use Control_Mod
  use Work_Mod,    only: phi_x       => r_cell_01,  &
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
  type(Grid_Type) :: grid
  real            :: dt
  integer         :: ini
  type(Var_Type)  :: phi
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter
  real              :: f_ex, f_im
  real              :: phis
  real              :: a0, a12, a21
  real              :: ini_res, tol
  real              :: vis_eff
  real              :: phix_f, phiy_f, phiz_f
  real              :: vis_t_f
  character(len=80) :: precond
  integer           :: adv_scheme    ! space discratization advection (scheme)
  real              :: blend         ! blending coeff (1.0 central; 0. upwind)
  integer           :: td_inertia    ! time-disretization for inerita  
  integer           :: td_advection  ! time-disretization for advection
  integer           :: td_diffusion  ! time-disretization for diffusion 
  integer           :: td_cross_diff ! time-disretization for cross-difusion
  real              :: urf           ! under-relaxation factor                 
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

  a % val = 0.

  b(:) = 0.

  !-------------------------------------! 
  !   Initialize variables and fluxes   !
  !-------------------------------------! 

  call Control_Mod_Time_Integration_For_Inertia(td_inertia)
  call Control_Mod_Time_Integration_For_Advection(td_advection)
  call Control_Mod_Time_Integration_For_Diffusion(td_diffusion)
  call Control_Mod_Time_Integration_For_Cross_Diffusion(td_cross_diff)

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      phi % oo(c)   = phi % o(c)
      phi % o (c)   = phi % n(c)
      phi % a_oo(c) = phi % a_o(c)
      phi % a_o (c) = 0. 
      phi % d_oo(c) = phi % d_o(c)
      phi % d_o (c) = 0. 
      phi % c_oo(c) = phi % c_o(c)
      phi % c_o (c) = phi % c(c) 
    end do
  end if

  ! Gradients
  call Grad_Mod_For_Phi(grid, phi % n, 1, phi_x, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 2, phi_y, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 3, phi_z, .true.)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Turbulence(adv_scheme)
  call Control_Mod_Blending_Coefficient_For_Turbulence(blend)

  ! Compute phimax and phimin
  if(adv_scheme .ne. CENTRAL) then
    call Calculate_Minimum_Maximum(grid, phi % n, phi_min, phi_max)
    goto 1
  end if

  ! New values
1 do c = 1, grid % n_cells
    phi % a(c)    = 0.
    phi % c(c)    = 0.
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
      if(adv_scheme .ne. CENTRAL) then
        call Advection_Scheme(grid, phis, s, phi % n, phi_min, phi_max,  &
                              phi_x, phi_y, phi_z,                       &
                              grid % dx, grid % dy, grid % dz,           &
                              adv_scheme, blend) 
      end if 

      ! Compute advection term
      if(ini .eq. 1) then 
        if(c2  > 0) then
          phi % a_o(c1) = phi % a_o(c1) - flux(s)*phis
          phi % a_o(c2) = phi % a_o(c2) + flux(s)*phis
        else
          phi % a_o(c1) = phi % a_o(c1) - flux(s)*phis
        end if 
      end if
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

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashforth scheeme for convective fluxes
  if(td_advection .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c) = b(c) + (1.5 * phi % a_o(c) - 0.5 * phi % a_oo(c) - phi % c(c))
    end do  
  end if

  ! Crank-Nicholson scheeme for convective fluxes
  if(td_advection .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c) = b(c) + (0.5 * ( phi % a(c) + phi % a_o(c) ) - phi % c(c))
    end do  
  end if

  ! Fully implicit treatment of convective fluxes 
  if(td_advection .eq. FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      b(c) = b(c) + (phi % a(c) - phi % c(c))
    end do  
  end if     
            
  ! New values
  do c = 1, grid % n_cells
    phi % c(c) = 0.
  end do
  
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
    vis_t_f = fw(s)*vis_t(c1) + (1.0-fw(s))*vis_t(c2)

    vis_eff = viscosity + vis_t_f

    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      if(turbulence_model_variant .ne. STABILIZED) then
        vis_eff = 1.5*viscosity 
      end if
    end if

    phix_f = fw(s) *phi_x(c1) + (1.0 - fw(s)) * phi_x(c2)
    phiy_f = fw(s) *phi_y(c1) + (1.0 - fw(s)) * phi_y(c2)
    phiz_f = fw(s) *phi_z(c1) + (1.0 - fw(s)) * phi_z(c2)


    ! Total (exact) diffusive flux plus turb. diffusion
    f_ex = vis_eff * (  phix_f * grid % sx(s)  &
                    + phiy_f * grid % sy(s)  &
                    + phiz_f * grid % sz(s) ) 

    a0 = vis_eff * f_coef(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: f_coef is
    !  not corrected at interface between materials)
    f_im=( phix_f*grid % dx(s)                      &
         +phiy_f*grid % dy(s)                      &
         +phiz_f*grid % dz(s))*a0

    ! Straight diffusion part 
    if(ini .eq. 1) then
      if(c2  > 0) then
        phi % d_o(c1) = phi % d_o(c1) + (phi % n(c2)-phi % n(c1))*a0   
        phi % d_o(c2) = phi % d_o(c2) - (phi % n(c2)-phi % n(c1))*a0    
      else
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. SYMMETRY) then
          phi % d_o(c1) = phi % d_o(c1) + (phi % n(c2)-phi % n(c1))*a0   
        end if 
      end if 
    end if

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex - f_im 
    if(c2  > 0) then
      phi % c(c2) = phi % c(c2) - f_ex + f_im 
    end if 

    ! Compute coefficients for the sysytem matrix
    if( (td_diffusion .eq. CRANK_NICOLSON) .or.  &
        (td_diffusion .eq. FULLY_IMPLICIT) ) then  

      if(td_diffusion .eq. CRANK_NICOLSON) then
        a12 = 0.5 * a0 
        a21 = 0.5 * a0 
      end if

      if(td_diffusion .eq. FULLY_IMPLICIT) then
        a12 = a0 
        a21 = a0 
      end if

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
!!!        (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT).or.    &
           (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) ) then
          a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
          b(c1) = b(c1) + a12 * phi % n(c2)
        end if
      end if     

    end if

  end do  ! through faces

  !------------------------------!
  !   Turbulent diffusion term   !
  !------------------------------!
  if(turbulence_model_variant .ne. STABILIZED) then
    if(phi % name .eq. 'EPS') then
      c_mu_d = 0.18        
    else
      c_mu_d = 0.22
    end if 
    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then        
      do c = 1, grid % n_cells
        u1uj_phij(c) = c_mu_d / phi % sigma * kin % n(c) / eps % n(c)  &
                     * (  uu % n(c) * phi_x(c)                         &
                        + uv % n(c) * phi_y(c)                         &
                        + uw % n(c) * phi_z(c))                        &
                     - viscosity * phi_x(c)

        u2uj_phij(c) = c_mu_d / phi % sigma * kin % n(c) / eps % n(c)  &
                     * (  uv % n(c) * phi_x(c)                         &
                        + vv % n(c) * phi_y(c)                         &
                        + vw % n(c) * phi_z(c))                        &
                     - viscosity * phi_y(c)

        u3uj_phij(c) = c_mu_d / phi % sigma * kin % n(c) / eps % n(c)  &
                     * (  uw % n(c) * phi_x(c)                         &
                        + vw % n(c) * phi_y(c)                         &
                        + ww % n(c) * phi_z(c))                        &
                     - viscosity * phi_z(c)
      end do
    else if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      do c = 1, grid % n_cells
        u1uj_phij(c) = c_mu_d / phi % sigma * t_scale(c)  &
                     * (  uu % n(c) * phi_x(c)            &
                        + uv % n(c) * phi_y(c)            &
                        + uw % n(c) * phi_z(c)) 

        u2uj_phij(c) = c_mu_d / phi % sigma * t_scale(c)  &
                     * (  uv % n(c) * phi_x(c)            &
                        + vv % n(c) * phi_y(c)            &
                        + vw % n(c) * phi_z(c)) 

        u3uj_phij(c) = c_mu_d / phi % sigma * t_scale(c)  &
                     * (  uw % n(c) * phi_x(c)            &
                        + vw % n(c) * phi_y(c)            &
                        + ww % n(c) * phi_z(c)) 
      end do
    end if
    call Grad_Mod_For_Phi(grid, u1uj_phij, 1, u1uj_phij_x, .true.)
    call Grad_Mod_For_Phi(grid, u2uj_phij, 2, u2uj_phij_y, .true.)
    call Grad_Mod_For_Phi(grid, u3uj_phij, 3, u3uj_phij_z, .true.)

    do c = 1, grid % n_cells
      b(c) = b(c) + (  u1uj_phij_x(c)  &
                     + u2uj_phij_y(c)  &
                     + u3uj_phij_z(c) ) * grid % vol(c)
    end do

    !------------------------------------------------------------------!
    !   Here we clean up transport equation from the false diffusion   !
    !------------------------------------------------------------------!
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .and.  &
       turbulence_model_variant .ne. STABILIZED) then
      do s = 1, grid % n_faces

        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        vis_eff = (fw(s)*vis_t(c1)+(1.0-fw(s))*vis_t(c2)) 

        phix_f = fw(s)*phi_x(c1) + (1.0-fw(s))*phi_x(c2)
        phiy_f = fw(s)*phi_y(c1) + (1.0-fw(s))*phi_y(c2)
        phiz_f = fw(s)*phi_z(c1) + (1.0-fw(s))*phi_z(c2)
        f_ex = vis_eff * (  phix_f * grid % sx(s)  &
                         + phiy_f * grid % sy(s)  &
                         + phiz_f * grid % sz(s))
        a0 = vis_eff * f_coef(s)
        f_im = (   phix_f * grid % dx(s)      &
                + phiy_f * grid % dy(s)      &
                + phiz_f * grid % dz(s)) * a0

        b(c1) = b(c1)                                            &
              - vis_eff * (phi % n(c2) - phi%n(c1)) * f_coef(s)  &
              - f_ex + f_im
        if(c2  > 0) then
          b(c2) = b(c2)                                            &
                + vis_eff * (phi % n(c2) - phi%n(c1)) * f_coef(s)  &
                + f_ex - f_im
        end if
      end do
    end if
  end if

  !---------------------------------!
  !     Temporal discretization     !
  !---------------------------------!

  ! Adams-Bashfort scheeme for diffusion fluxes
  if(td_diffusion .eq. ADAMS_BASHFORTH) then 
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * phi % d_o(c) - 0.5 * phi % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(td_diffusion .eq. CRANK_NICOLSON) then 
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * phi % d_o(c)
    end do  
  end if
                 
  ! Fully implicit treatment for difusive terms
  !  is handled via the linear system of equations 

  ! Adams-Bashfort scheeme for cross diffusion 
  if(td_cross_diff .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * phi % c_o(c) - 0.5 * phi % c_oo(c)
    end do 
  end if
    
  ! Crank-Nicholson scheme for cross difusive terms
  if(td_cross_diff .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * phi % c(c) + 0.5 * phi % c_o(c)
    end do 
  end if

  ! Fully implicit treatment for cross difusive terms
  if(td_cross_diff .eq. FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      b(c) = b(c) + phi % c(c) 
    end do 
  end if

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  ! Two time levels; Linear interpolation
  if(td_inertia .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = density*grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c) = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_inertia .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = density*grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0*a0 * phi % o(c) - 0.5*a0 * phi % oo(c)
    end do
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then 
    call Grad_Mod_For_Phi(grid, f22 % n, 1, f22 % x, .true.) ! df22/dx
    call Grad_Mod_For_Phi(grid, f22 % n, 2, f22 % y, .true.) ! df22/dy
    call Grad_Mod_For_Phi(grid, f22 % n, 3, f22 % z, .true.) ! df22/dz

    call Source_EBM              (grid, phi % name)
  else
    call Source_Hanjalic_Jakirlic(grid, phi % name)        
  end if                

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !    
  !---------------------------------!

  ! Set under-relaxation factor
  urf = 1.0
  call Control_Mod_Simple_Underrelaxation_For_Turbulence(urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - urf)*phi % n(c) / urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / urf
!?????? Asave(c) = a % val(a % dia(c)) ??????
  end do

  call Control_Mod_Tolerance_For_Turbulence_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set the default value for number of iterations
  niter = 6

  ! Over-ride if specified in control file
  call Control_Mod_Max_Iterations_For_Turbulence_Solver(niter)

  call Cg(a, phi % n, b, precond, niter, tol, ini_res, phi % res)

  if( phi % name .eq. 'UU' )   &
    call Info_Mod_Iter_Fill_At(3, 1, phi % name, niter, phi % res)
  if( phi % name .eq. 'VV' )   &
    call Info_Mod_Iter_Fill_At(3, 2, phi % name, niter, phi % res)
  if( phi % name .eq. 'WW' )   &
    call Info_Mod_Iter_Fill_At(3, 3, phi % name, niter, phi % res)
  if( phi % name .eq. 'UV' )   &
    call Info_Mod_Iter_Fill_At(4, 1, phi % name, niter, phi % res)
  if( phi % name .eq. 'UW' )   &
    call Info_Mod_Iter_Fill_At(4, 2, phi % name, niter, phi % res)
  if( phi % name .eq. 'VW' )   &
    call Info_Mod_Iter_Fill_At(4, 3, phi % name, niter, phi % res)
  if( phi % name .eq. 'EPS' )  &
    call Info_Mod_Iter_Fill_At(4, 4, phi % name, niter, phi % res)

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
