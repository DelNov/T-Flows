!==============================================================================!
  subroutine Compute_Energy(grid, dt, ini, phi)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
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
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Grid_Type) :: grid
  integer         :: ini
  real            :: dt
  type(Var_Type)  :: phi
!----------------------------------[Calling]-----------------------------------!
  real :: Turbulent_Prandtl_Number
!-----------------------------------[Locals]-----------------------------------! 
  integer           :: n, c, s, c1, c2, niter
  real              :: a0, a12, a21
  real              :: ini_res, tol
  real              :: con_eff1, f_ex1, f_im1, phix_f1, phiy_f1, phiz_f1
  real              :: con_eff2, f_ex2, f_im2, phix_f2, phiy_f2, phiz_f2
  real              :: phis, pr_t1, pr_t2
  real              :: ut_s, vt_s, wt_s, t_stress, con_t
  character(len=80) :: precond       ! preconditioner
  integer           :: adv_scheme    ! space-discretiztion of advection scheme)
  real              :: blend         ! blending coeff (1.0 central; 0.0 upwind)
  integer           :: td_inertia    ! time-disretization for inerita  
  integer           :: td_advection  ! time-disretization for advection
  integer           :: td_diffusion  ! time-disretization for diffusion 
  integer           :: td_cross_diff ! time-disretization for cross-difusion
  real              :: urf           ! under-relaxation factor
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

  do n = 1, a % row(grid % n_cells+1)  ! this is number of non-zeros plus 1
    a % val(n) = 0.0
  end do
  a % val = 0.0

  b(:) = 0.0


  !-------------------------------------! 
  !   Initialize variables and fluxes   !
  !-------------------------------------! 

  call Control_Mod_Time_Integration_For_Inertia(td_inertia)
  call Control_Mod_Time_Integration_For_Advection(td_advection)
  call Control_Mod_Time_Integration_For_Diffusion(td_diffusion)
  call Control_Mod_Time_Integration_For_Cross_Diffusion(td_cross_diff)

  ! Old values (o and oo)
  if(ini.lt.2) then
    do c = 1, grid % n_cells
      phi % oo(c)   = phi % o(c)
      phi % o (c)   = phi % n(c)
      phi % a_oo(c) = phi % a_o(c)
      phi % a_o (c) = 0.0
      phi % d_oo(c) = phi % d_o(c)
      phi % d_o (c) = 0.0 
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
  call Control_Mod_Advection_Scheme_For_Energy(adv_scheme)
  call Control_Mod_Blending_Coefficient_For_Energy(blend)

  ! Compute phimax and phimin
  if(adv_scheme .ne. CENTRAL) then
    call Calculate_Minimum_Maximum(grid, phi % n, phi_min, phi_max)
    goto 1  ! why this???
  end if

  ! New values
1 do c = 1, grid % n_cells
    phi % a(c) = 0.0
    phi % c(c) = 0.0  ! use phi % c for upwind advective fluxes
  end do

  !----------------------------------!
  !   Browse through all the faces   !
  !----------------------------------!
  do s=1,grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    phis =      grid % f(s)  * phi % n(c1)   &
         + (1.0-grid % f(s)) * phi % n(c2)

    ! Compute phis with desired advection scheme
    if(adv_scheme .ne. CENTRAL) then
      call Advection_Scheme(grid, phis, s, phi % n, phi_min, phi_max,  &
                            phi_x, phi_y, phi_z,                       &
                            grid % dx, grid % dy, grid % dz,           &
                            adv_scheme, blend) 
    end if

    ! Compute advection term
    if(ini.eq.1) then
      if(c2.gt.0) then
        phi % a_o(c1) = phi % a_o(c1)-flux(s)*phis*capacity
        phi % a_o(c2) = phi % a_o(c2)+flux(s)*phis*capacity
      else
        phi % a_o(c1) = phi % a_o(c1)-flux(s)*phis*capacity
      end if
    end if
    if(c2.gt.0) then
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
      phi % a(c2) = phi % a(c2)+flux(s)*phis*capacity
    else
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
    end if

    ! Store upwinded part of the advection term in "c"
    if(pressure_momentum_coupling .ne. PROJECTION) then
      if(flux(s).lt.0) then   ! from c2 to c1
        phi % c(c1) = phi % c(c1)-flux(s)*phi % n(c2) * capacity
        if(c2.gt.0) then
          phi % c(c2) = phi % c(c2)+flux(s)*phi % n(c2) * capacity
        end if
      else
        phi % c(c1) = phi % c(c1)-flux(s)*phi % n(c1) * capacity
        if(c2.gt.0) then
          phi % c(c2) = phi % c(c2)+flux(s)*phi % n(c1) * capacity
        end if
      end if
    end if  

 
  end do  ! through sides

  !-----------------------------!
  !                             !
  !   Temporal discretization   !
  !                             !
  !-----------------------------!
   
  ! Adams-Bashforth scheeme for advection fluxes
  if(td_advection .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5*phi % a_o(c) - 0.5*phi % a_oo(c) - phi % c(c)
    end do
  end if

  ! Crank-Nicholson scheeme for advection fluxes
  if(td_advection .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * (phi % a(c) + phi % a_o(c)) - phi % c(c)
    end do
  end if

  ! Fully implicit treatment of advection fluxes
  if(td_advection .eq. FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      b(c) = b(c) + phi % a(c) - phi % c(c)
    end do
  end if

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  ! Set phi % c back to zero 
  do c = 1, grid % n_cells
    phi % c(c) = 0.0  
  end do

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  if(turbulence_model .ne. NONE .and.  &
     turbulence_model .ne. DNS) then
    call Control_Mod_Turbulent_Prandtl_Number(pr_t)  ! get default pr_t (0.9)
  end if

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(turbulence_model .ne. LES_SMAGORINSKY .and.  &
       turbulence_model .ne. LES_DYNAMIC     .and.  &
       turbulence_model .ne. LES_WALE        .and.  &
       turbulence_model .ne. NONE            .and.  &
       turbulence_model .ne. DNS) then
      pr_t1 = Turbulent_Prandtl_Number(grid, c1)
      pr_t2 = Turbulent_Prandtl_Number(grid, c2)
      pr_t  = fw(s) * pr_t1 + (1.0 - fw(s)) * pr_t2
    end if

    ! Gradients on the cell face (fw corrects situation close to the wall)
    phix_f1 = fw(s)*phi_x(c1) + (1.0-fw(s))*phi_x(c2) 
    phiy_f1 = fw(s)*phi_y(c1) + (1.0-fw(s))*phi_y(c2)
    phiz_f1 = fw(s)*phi_z(c1) + (1.0-fw(s))*phi_z(c2)
    phix_f2 = phix_f1
    phiy_f2 = phiy_f1
    phiz_f2 = phiz_f1
    if(turbulence_model .ne. NONE .and.  &
       turbulence_model .ne. DNS) then
      con_eff1 =        fw(s)  * (conductivity+capacity*vis_t(c1)/pr_t)  &
               + (1.0 - fw(s)) * (conductivity+capacity*vis_t(c2)/pr_t)
      con_t    = fw(s) * capacity*vis_t(c1)/pr_t &
               + (1.0 - fw(s)) * capacity*vis_t(c2)/pr_t
    else
      con_eff1 = conductivity
    end if

    con_eff2 = con_eff1 

    if(turbulence_model .eq. K_EPS        .or.  &
       turbulence_model .eq. K_EPS_ZETA_F .or.  &
       turbulence_model .eq. HYBRID_LES_RANS) then
      if(c2 < 0) then
        if(Var_Mod_Bnd_Cell_Type(phi, c2) .eq. WALL .or.  &
           Var_Mod_Bnd_Cell_Type(phi, c2) .eq. WALLFL) then
          con_eff1 = con_wall(c1)
          con_eff2 = con_eff1
        end if
      end if
    end if

    ! Total (exact) diffusive flux
    f_ex1 = con_eff1 * (  phix_f1 * grid % sx(s)   &
                        + phiy_f1 * grid % sy(s)   &
                        + phiz_f1 * grid % sz(s))
    f_ex2 = con_eff2 * (  phix_f2 * grid % sx(s)   &
                        + phiy_f2 * grid % sy(s)   &
                        + phiz_f2 * grid % sz(s))

    ! Implicit diffusive flux
    f_im1 = con_eff1 * f_coef(s)           &
          * (  phix_f1 * grid % dx(s)      &
             + phiy_f1 * grid % dy(s)      &
             + phiz_f1 * grid % dz(s) )
    f_im2 = con_eff2 * f_coef(s)           &
          * (  phix_f2 * grid % dx(s)      &
             + phiy_f2 * grid % dy(s)      &
             + phiz_f2 * grid % dz(s) )

    ! Straight diffusion part 
    if(ini.lt.2) then
      if(c2.gt.0) then
        phi % d_o(c1) = phi % d_o(c1)  &
                      + con_eff1*f_coef(s)*(phi % n(c2) - phi % n(c1))
        phi % d_o(c2) = phi % d_o(c2)  &
                      - con_eff2*f_coef(s)*(phi % n(c2) - phi % n(c1))
      else
        if(Var_Mod_Bnd_Cell_Type(phi, c2) .ne. SYMMETRY) then 
          phi % d_o(c1) = phi % d_o(c1)  &
                        + con_eff1*f_coef(s)*(phi % n(c2) - phi % n(c1))
        end if
      end if
    end if

    ! Turbulent heat fluxes according to GGDH scheme
    ! (first line is GGDH, second line is SGDH substratced 
    if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
       turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
      ut_s =  (     grid % f(s)  * ut % n(c1)  &
           +  (1. - grid % f(s)) * ut % n(c2))
      vt_s =  (     grid % f(s)  * vt % n(c1)  &
           +  (1. - grid % f(s)) * vt % n(c2))
      wt_s =  (     grid % f(s)  * wt % n(c1)  &
           +  (1. - grid % f(s)) * wt % n(c2))
      t_stress = - (  ut_s * grid % sx(s)                    &
                    + vt_s * grid % sy(s)                    &
                    + wt_s * grid % sz(s) )                  &
                    - (con_t * (  phix_f1 * grid % sx(s)     &
                                + phiy_f1 * grid % sy(s)     &
                                + phiz_f1 * grid % sz(s)) )
    end if

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex1 - f_im1
    if(c2.gt.0) then
      phi % c(c2) = phi % c(c2) - f_ex2 + f_im2
    end if

    ! Calculate the coefficients for the sysytem matrix
    if( (td_diffusion .eq. CRANK_NICOLSON) .or.  &
        (td_diffusion .eq. FULLY_IMPLICIT) ) then

      if(td_diffusion .eq. CRANK_NICOLSON) then
        a12 = .5 * con_eff1 * f_coef(s)
        a21 = .5 * con_eff2 * f_coef(s)
      end if

      if(td_diffusion .eq. FULLY_IMPLICIT) then 
        a12 = con_eff1 * f_coef(s)
        a21 = con_eff2 * f_coef(s)
      end if

      if(pressure_momentum_coupling .ne. PROJECTION) then
        a12 = a12  - min(flux(s), 0.0) * capacity
        a21 = a21  + max(flux(s), 0.0) * capacity
      end if

      ! Fill the system matrix
      if(c2.gt.0) then
        a % val(a % dia(c1))  = a % val(a % dia(c1)) + a12
        a % val(a % dia(c2))  = a % val(a % dia(c2)) + a21
        a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
        a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
      else if(c2.lt.0) then
        ! Outflow is included because of the flux 
        ! corrections which also affects velocities
        if( (Var_Mod_Bnd_Cell_Type(phi, c2) .eq. INFLOW) .or.  &
            (Var_Mod_Bnd_Cell_Type(phi, c2) .eq. WALL)   .or.  &
            (Var_Mod_Bnd_Cell_Type(phi, c2) .eq. CONVECT) ) then
          a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
          b(c1)  = b(c1)  + a12 * phi % n(c2)
        ! In case of wallflux 
        else if(Var_Mod_Bnd_Cell_Type(phi, c2) .eq. WALLFL) then
          b(c1) = b(c1) + grid % s(s) * phi % q(c2)
        end if 
      end if

    end if

  end do  ! through sides

  !-----------------------------!
  !                             !
  !   Temporal discretization   !
  !                             !
  !-----------------------------!
   
  ! Adams-Bashfort scheeme for diffusion fluxes
  if(td_diffusion .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c)  = b(c) + 1.5*phi % d_o(c) - 0.5*phi % d_oo(c)
    end do  
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(td_diffusion .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c)  = b(c) + 0.5*phi % d_o(c)
    end do  
  end if

  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations

  ! Adams-Bashfort scheeme for cross diffusion 
  if(td_cross_diff .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c)  = b(c) + 1.5*phi % c_o(c) - 0.5*phi % c_oo(c)
    end do 
  end if

  ! Crank-Nicholson scheme for cross difusive terms
  if(td_cross_diff .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      if( (phi % c(c)+phi % c_o(c))  >= 0) then
        b(c)  = b(c) + 0.5*(phi % c(c) + phi % c_o(c))
      else
        a % val(a % dia(c)) = a % val(a % dia(c)) &
             - 0.5 * (phi % c(c) + phi % c_o(c)) / (phi % n(c) + MICRO)
      end if
    end do
  end if
 
  ! Fully implicit treatment for cross difusive terms
  if(td_cross_diff .eq. FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      if(phi % c(c) >= 0) then
        b(c)  = b(c) + phi % c(c)
      else
        a % val(a % dia(c)) = a % val(a % dia(c))  &
                            - phi % c(c) / (phi % n(c) + MICRO)
      end if
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
      a0 = capacity * density * grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c)  = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_inertia .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = capacity * density * grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c)  = b(c) + 2.0 * a0 * phi % o(c) - 0.5 * a0 * phi % oo(c)
    end do
  end if

  call User_Mod_Source(grid, phi, a, b)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Type of coupling is important
  call Control_Mod_Pressure_Momentum_Coupling()

  ! Set under-relaxation factor
  urf = 1.0
  if(pressure_momentum_coupling .eq. SIMPLE)   &
    call Control_Mod_Simple_Underrelaxation_For_Energy(urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - urf) * phi % n(c)  &
                                      / urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / urf
  end do  

  ! Get solver tolerance
  call Control_Mod_Tolerance_For_Energy_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set number of iterations based on coupling scheme
  if(pressure_momentum_coupling .eq. PROJECTION) niter = 10
  if(pressure_momentum_coupling .eq. SIMPLE)     niter =  5

  ! Over-ride if specified in control file
  call Control_Mod_Max_Iterations_For_Energy_Solver(niter)

  call Bicg(a,        &
            phi % n,  &
            b,        &
            precond,  &
            niter,    &
            tol,      &
            ini_res,  &
            phi % res)

  call Info_Mod_Iter_Fill_At(2, 4, phi % name, niter, phi % res)

  call Comm_Mod_Exchange_Real(grid, phi % n)

  end subroutine
