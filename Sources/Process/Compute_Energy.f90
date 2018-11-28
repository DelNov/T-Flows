!==============================================================================!
  subroutine Compute_Energy(grid, sol, dt, ini, phi)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Comm_Mod
  use Var_Mod,      only: Var_Type
  use Grid_Mod,     only: Grid_Type
  use Grad_Mod
  use Info_Mod
  use Numerics_Mod, only: CENTRAL, LINEAR, PARABOLIC
  use Solver_Mod,   only: Solver_Type, Bicg, Cg, Cgs
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
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Grid_Type)           :: grid
  type(Solver_Type), target :: sol
  integer                   :: ini
  real                      :: dt
  type(Var_Type)            :: phi
!----------------------------------[Calling]-----------------------------------!
  real :: Turbulent_Prandtl_Number
!-----------------------------------[Locals]-----------------------------------! 
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
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
  integer           :: td_scheme     ! time-disretization for inerita  
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

  ! Take aliases
  a => sol % a
  b => sol % b

  do n = 1, a % row(grid % n_cells+1)  ! this is number of non-zeros plus 1
    a % val(n) = 0.0
  end do
  a % val = 0.0

  b(:) = 0.0

  ! Old values (o and oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
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
  end if

  ! New values
  do c = 1, grid % n_cells
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
    if(c2 > 0) then
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
      phi % a(c2) = phi % a(c2)+flux(s)*phis*capacity
    else
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
    end if

    ! Store upwinded part of the advection term in "c"
    if(flux(s).lt.0) then   ! from c2 to c1
      phi % c(c1) = phi % c(c1)-flux(s)*phi % n(c2) * capacity
      if(c2 > 0) then
        phi % c(c2) = phi % c(c2)+flux(s)*phi % n(c2) * capacity
      end if
    else
      phi % c(c1) = phi % c(c1)-flux(s)*phi % n(c1) * capacity
      if(c2 > 0) then
        phi % c(c2) = phi % c(c2)+flux(s)*phi % n(c1) * capacity
      end if
    end if

  end do  ! through faces

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + phi % a(c) - phi % c(c)
  end do

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

    ! Total (exact) diffusion flux
    f_ex1 = con_eff1 * (  phix_f1 * grid % sx(s)   &
                        + phiy_f1 * grid % sy(s)   &
                        + phiz_f1 * grid % sz(s))
    f_ex2 = con_eff2 * (  phix_f2 * grid % sx(s)   &
                        + phiy_f2 * grid % sy(s)   &
                        + phiz_f2 * grid % sz(s))

    ! Implicit diffusion flux
    f_im1 = con_eff1 * f_coef(s)           &
          * (  phix_f1 * grid % dx(s)      &
             + phiy_f1 * grid % dy(s)      &
             + phiz_f1 * grid % dz(s) )
    f_im2 = con_eff2 * f_coef(s)           &
          * (  phix_f2 * grid % dx(s)      &
             + phiy_f2 * grid % dy(s)      &
             + phiz_f2 * grid % dz(s) )

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex1 - f_im1 
    if(c2 > 0) then
      phi % c(c2) = phi % c(c2) - f_ex2 + f_im2 
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
      ut_s =  (     fw(s)  * ut % n(c1)  &
           +  (1. - fw(s)) * ut % n(c2))
      vt_s =  (     fw(s)  * vt % n(c1)  &
           +  (1. - fw(s)) * vt % n(c2))
      wt_s =  (     fw(s)  * wt % n(c1)  &
           +  (1. - fw(s)) * wt % n(c2))
      t_stress = - (  ut_s * grid % sx(s)                    &
                    + vt_s * grid % sy(s)                    &
                    + wt_s * grid % sz(s) )                  &
                    - (con_t * (  phix_f1 * grid % sx(s)     &
                                + phiy_f1 * grid % sy(s)     &
                                + phiz_f1 * grid % sz(s)) )

      ! Put the influence of turbulent heat fluxes explicitly in the system
      b(c1) = b(c1) + t_stress
      if(c2 > 0) then
        b(c2) = b(c2) - t_stress
      end if
    end if  ! if models are of RSM type

    ! Calculate the coefficients for the sysytem matrix
    a12 = con_eff1 * f_coef(s)
    a21 = con_eff2 * f_coef(s)

    a12 = a12  - min(flux(s), 0.0) * capacity
    a21 = a21  + max(flux(s), 0.0) * capacity

    ! Fill the system matrix
    if(c2 > 0) then
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

  end do  ! through sides

  ! Cross diffusion terms are treated explicity
  do c = 1, grid % n_cells
    if(phi % c(c) >= 0) then
      b(c)  = b(c) + phi % c(c)
    else
      a % val(a % dia(c)) = a % val(a % dia(c))  &
                          - phi % c(c) / (phi % n(c) + MICRO)
    end if
  end do


  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  call Control_Mod_Time_Integration_Scheme(td_scheme)

  ! Two time levels; Linear interpolation
  if(td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = capacity * density * grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c)  = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_scheme .eq. PARABOLIC) then
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

  ! Set under-relaxation factor, then overwrite with control file if specified
  urf = 0.7
  call Control_Mod_Simple_Underrelaxation_For_Energy(urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - urf) * phi % n(c) / urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / urf
  end do

  ! Get solver tolerance
  call Control_Mod_Tolerance_For_Energy_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set number of iterations then overwrite with control file if specified
  niter =  5
  call Control_Mod_Max_Iterations_For_Energy_Solver(niter)

  call Bicg(sol,      &
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
