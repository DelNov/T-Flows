!==============================================================================!
  subroutine User_Mod_Compute_Scalar(grid, sol, dt, ini, phi)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for use scalar.                          !
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
!-----------------------------------[Arguments]--------------------------------!
  type(Grid_Type)   :: grid
  type(Solver_Type), target :: sol
  integer         :: ini
  real            :: dt
  type(Var_Type)  :: phi
!----------------------------------[Calling]-----------------------------------!
  real :: Turbulent_Prandtl_Number
!-----------------------------------[Locals]-----------------------------------! 
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer           :: n, c, s, c1, c2, niter, mat, row, col
  real              :: a0, a12, a21
  real              :: ini_res, tol, ns
  real              :: con_eff1, f_ex1, f_im1, phix_f1, phiy_f1, phiz_f1
  real              :: con_eff2, f_ex2, f_im2, phix_f2, phiy_f2, phiz_f2
  real              :: phis, pr_t1, pr_t2
  character(len=80) :: precond    ! preconditioner
  integer           :: adv_scheme  ! space-discretiztion of advection scheme)
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

  ! Take aliases
  a => sol % a
  b => sol % b % val

  do n = 1, a % row(grid % n_cells+1) ! to je broj nonzero + 1
    a % val(n) = 0.0
  end do
  a % val = 0.0

  b(:) = 0.0

  !-------------------------------------!
  !   Initialize variables and fluxes   !
  !-------------------------------------!

  ! Old values (o and oo)
  if(ini.lt.2) then
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
  call Control_Mod_Advection_Scheme_For_User_Scalars(adv_scheme)
  call Control_Mod_Blending_Coefficient_For_User_Scalars(blend)

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
    if(c2 > 0) then
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
      phi % a(c2) = phi % a(c2)+flux(s)*phis*capacity
    else
      phi % a(c1) = phi % a(c1)-flux(s)*phis*capacity
    end if

    ! Store upwinded part of the advection term in "c"
    if(flux(s) < 0.0) then   ! from c2 to c1
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

  end do  ! through sides

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!
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
  call Control_Mod_Turbulent_Prandtl_Number(pr_t)  ! get default pr_t (0.9)

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(turbulence_model .ne. LES_SMAGORINSKY .or.  &
       turbulence_model .ne. LES_DYNAMIC     .or.  &
       turbulence_model .ne. LES_WALE        .or.  &
       turbulence_model .ne. DNS) then
      pr_t1 = Turbulent_Prandtl_Number(grid, c1)
      pr_t2 = Turbulent_Prandtl_Number(grid, c2)
      pr_t  = grid % fw(s) * pr_t1 + (1.0-grid % fw(s)) * pr_t2
    end if

    ! Gradients on the cell face 
    if(c2 > 0) then
      phix_f1 = grid % fw(s)*phi_x(c1) + (1.0-grid % fw(s))*phi_x(c2)
      phiy_f1 = grid % fw(s)*phi_y(c1) + (1.0-grid % fw(s))*phi_y(c2)
      phiz_f1 = grid % fw(s)*phi_z(c1) + (1.0-grid % fw(s))*phi_z(c2)
      phix_f2 = phix_f1 
      phiy_f2 = phiy_f1 
      phiz_f2 = phiz_f1 
      con_eff1 =     grid % f(s)  * (conductivity+capacity*vis_t(c1)/pr_t)  &
               + (1.-grid % f(s)) * (conductivity+capacity*vis_t(c2)/pr_t)
      con_eff2 = con_eff1 
    else
      phix_f1 = phi_x(c1) 
      phiy_f1 = phi_y(c1) 
      phiz_f1 = phi_z(c1) 
      phix_f2 = phix_f1 
      phiy_f2 = phiy_f1 
      phiz_f2 = phiz_f1 
      con_eff1 = conductivity + capacity * vis_t(c1) / pr_t   
      con_eff2 = con_eff1 
    end if


    if(turbulence_model .eq. K_EPS .or.  &
       turbulence_model .eq. K_EPS_ZETA_F) then 
      if(c2 < 0) then
        if(Var_Mod_Bnd_Cell_Type(phi,c2) .eq. WALL .or.  &
           Var_Mod_Bnd_Cell_Type(phi,c2) .eq. WALLFL) then
          con_eff1 = con_wall(c1)
          con_eff2 = con_eff1
        end if
      end if
    end if  

    ! Total (exact) diffusive flux
    f_ex1 = con_eff1 * (  phix_f1 * grid % sx(s)  &
                        + phiy_f1 * grid % sy(s)  &
                        + phiz_f1 * grid % sz(s))
    f_ex2 = con_eff2 * (  phix_f2 * grid % sx(s)  &
                        + phiy_f2 * grid % sy(s)  &
                        + phiz_f2 * grid % sz(s))

    ! Implicit diffusive flux
    f_im1 = con_eff1 * a % fc(s)          &
          * (  phix_f1 * grid % dx(s)      &
             + phiy_f1 * grid % dy(s)      &
             + phiz_f1 * grid % dz(s) )
    f_im2 = con_eff2 * a % fc(s)          &
          * (  phix_f2 * grid % dx(s)      &
             + phiy_f2 * grid % dy(s)      &
             + phiz_f2 * grid % dz(s) )

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex1 - f_im1 
    if(c2 .gt. 0) then
      phi % c(c2) = phi % c(c2) - f_ex2 + f_im2 
    end if 

    ! Calculate the coefficients for the sysytem matrix

    a12 = con_eff1 * a % fc(s)
    a21 = con_eff2 * a % fc(s)

    a12 = a12  - min(flux(s), 0.0) * capacity
    a21 = a21  + max(flux(s), 0.0) * capacity

    ! Fill the system matrix
    if(c2 > 0) then
      a % val(a % dia(c1))  = a % val(a % dia(c1)) + a12
      a % val(a % dia(c2))  = a % val(a % dia(c2)) + a21
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
    else if(c2 < 0) then

      ! Outflow is included because of the flux 
      ! corrections which also affects velocities
      if( (Var_Mod_Bnd_Cell_Type(phi,c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cell_Type(phi,c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cell_Type(phi,c2) .eq. CONVECT) ) then    
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * phi % n(c2)

      ! In case of wallflux 
      else if(Var_Mod_Bnd_Cell_Type(phi,c2) .eq. WALLFL) then
        b(c1) = b(c1) + grid % s(s) * phi % q(c2)
      end if

    end if

  end do  ! through sides

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations 

  ! Fully implicit treatment for cross difusive terms
  do c = 1, grid % n_cells
    if(phi % c(c) >= 0) then
      b(c)  = b(c) + phi % c(c)
    else
      a % val(a % dia(c)) = a % val(a % dia(c))  &
                          - phi % c(c)/(phi % n(c)+1.e-6)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  call Control_Mod_Time_Integration_Scheme(td_inertia)

  ! Two time levels; Linear interpolation
  if(td_inertia .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = capacity * density * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c)  = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_inertia .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = capacity * density * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c)  = b(c) + 2.0 * a0 * phi % o(c) - 0.5 * a0 * phi % oo(c)
    end do
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(turbulence_model_variant .ne. STABILIZED) then
      do c = 1, grid % n_cells
        u1uj_phij(c) = -0.22*t_scale(c) *&
                   (uu%n(c)*phi_x(c)+uv%n(c)*phi_y(c)+uw%n(c)*phi_z(c))
        u2uj_phij(c) = -0.22*t_scale(c)*&
                   (uv%n(c)*phi_x(c)+vv%n(c)*phi_y(c)+vw%n(c)*phi_z(c))
        u3uj_phij(c) = -0.22*t_scale(c)*&
                   (uw%n(c)*phi_x(c)+vw%n(c)*phi_y(c)+ww%n(c)*phi_z(c))
      end do
      call Grad_Mod_For_Phi(grid, u1uj_phij, 1, u1uj_phij_x, .true.)
      call Grad_Mod_For_Phi(grid, u2uj_phij, 2, u2uj_phij_y, .true.)
      call Grad_Mod_For_Phi(grid, u3uj_phij, 3, u3uj_phij_z, .true.)
      do c = 1, grid % n_cells
        b(c) = b(c) - (  u1uj_phij_x(c)  &
                       + u2uj_phij_y(c)  &
                       + u3uj_phij_z(c) ) * grid % vol(c)
      end do

      !------------------------------------------------------------------!
      !   Here we clean up transport equation from the false diffusion   !
      !------------------------------------------------------------------!
      do s = 1, grid % n_faces

        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        pr_t1 = Turbulent_Prandtl_Number(grid, c1)
        pr_t2 = Turbulent_Prandtl_Number(grid, c2)
        pr_t  = grid % fw(s) * pr_t1 + (1.0-grid % fw(s)) * pr_t2

        if(c2 > 0) then
          phix_f1 = grid % fw(s)*phi_x(c1) + (1.0-grid % fw(s))*phi_x(c2)
          phiy_f1 = grid % fw(s)*phi_y(c1) + (1.0-grid % fw(s))*phi_y(c2)
          phiz_f1 = grid % fw(s)*phi_z(c1) + (1.0-grid % fw(s))*phi_z(c2)
          phix_f2 = phix_f1 
          phiy_f2 = phiy_f1 
          phiz_f2 = phiz_f1 
          con_eff1 =      grid % f(s)  * (capacity*vis_t(c1)/pr_t )  &
                  + (1. - grid % f(s)) * (capacity*vis_t(c2)/pr_t )
          con_eff2 = con_eff1 
        else
          phix_f1 = phi_x(c1)
          phiy_f1 = phi_y(c1)
          phiz_f1 = phi_z(c1)
          phix_f2 = phix_f1 
          phiy_f2 = phiy_f1 
          phiz_f2 = phiz_f1 
          con_eff1 = capacity * vis_t(c1) / pr_t
          con_eff2 = con_eff1 
        end if

        ! Total (exact) diffusive flux
        f_ex1 = con_eff1 * (  phix_f1 * grid % sx(s)  &
                            + phiy_f1 * grid % sy(s)  &
                            + phiz_f1 * grid % sz(s))
        f_ex2 = con_eff2 * (  phix_f2 * grid % sx(s)  &
                            + phiy_f2 * grid % sy(s)  &
                            + phiz_f2 * grid % sz(s))

        ! Implicit diffusive flux
        f_im1 = con_eff1 * a % fc(s) *         &
                (  phix_f1 * grid % dx(s)      &
                 + phiy_f1 * grid % dy(s)      &
                 + phiz_f1 * grid % dz(s) )
        f_im2 = con_eff2 * a % fc(s) *         &
                (  phix_f2 * grid % dx(s)      &
                 + phiy_f2 * grid % dy(s)      &
                 + phiz_f2 * grid % dz(s) )

        b(c1) = b(c1) - con_eff1 * (phi % n(c2) - phi % n(c1)) * a % fc(s)  &
              - f_ex1 + f_im1
        if(c2  > 0) then
          b(c2) = b(c2) + con_eff1 * (phi % n(c2) - phi % n(c1)) * a % fc(s)  &
                + f_ex2 - f_im2
        end if
      end do
    end if  
  end if  

  call User_Mod_Source(grid, phi, a, b)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !    
  !---------------------------------!

  ! Set under-relaxation factor
  urf = 0.7
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

  ! Set number of iterations
  niter =  5
  call Control_Mod_Max_Iterations_For_Energy_Solver(niter)

  call Bicg(sol, phi % n, b, precond, niter, tol, ini_res, phi % res)

  read(phi % name(3:4), *) ns  ! reterive the number of scalar 
  row = ceiling(ns/4)          ! will be 1 (scal. 1-4), 2 (scal. 5-8), etc.
  col = ns - (row-1)*4         ! will be in range 1 - 4

  call Info_Mod_Iter_Fill_User_At(row, col, phi % name, niter, phi % res)
 
  call Comm_Mod_Exchange_Real(grid, phi % n)

  end subroutine
