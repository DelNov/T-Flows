!==============================================================================!
  subroutine Compute_Momentum(grid, sol, dt, ini, ui,  &
                              ui_i, ui_j, ui_k,   &
                              si, sj, sk,         &
                              di, dj, dk,         &
                              h_i, uj_i, uk_i)
!------------------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Comm_Mod
  use Var_Mod,      only: Var_Type
  use Grid_Mod,     only: Grid_Type
  use Bulk_Mod
  use Info_Mod,     only: Info_Mod_Iter_Fill_At
  use Numerics_Mod, only: CENTRAL, LINEAR, PARABOLIC
  use Solver_Mod,   only: Solver_Type, Bicg, Cg, Cgs
  use Matrix_Mod,   only: Matrix_Type
  use Control_Mod
  use User_Mod
  use Work_Mod,      only: ui_min  => r_cell_01,  &
                           ui_max  => r_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  type(Solver_Type), target :: sol
  real            :: dt
  integer         :: ini
  type(Var_Type)  :: ui
  real            :: ui_i(-grid % n_bnd_cells:grid % n_cells),  &
                     ui_j(-grid % n_bnd_cells:grid % n_cells),  &
                     ui_k(-grid % n_bnd_cells:grid % n_cells)
  real            :: si(grid % n_faces),  &
                     sj(grid % n_faces),  &
                     sk(grid % n_faces)
  real            :: di(grid % n_faces),  &
                     dj(grid % n_faces),  &
                     dk(grid % n_faces)
  real            :: h_i (-grid % n_bnd_cells:grid % n_cells),  &
                     uj_i(-grid % n_bnd_cells:grid % n_cells),  &
                     uk_i(-grid % n_bnd_cells:grid % n_cells)
  real            :: uu_f, vv_f, ww_f, uv_f, uw_f, vw_f
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer           :: s, c, c1, c2, niter
  real              :: f_ex, f_im, f_stress
  real              :: uis, vel_max
  real              :: a0, a12, a21
  real              :: ini_res, tol
  real              :: vis_eff, vis_tS
  real              :: ui_i_f,ui_j_f,ui_k_f,uj_i_f,uk_i_f
  character(len=80) :: precond
  integer           :: adv_scheme    ! space disretization of advection (scheme)
  real              :: blend         ! blending coeff (1.0 central; 0.0 upwind)
  integer           :: td_scheme     ! time-disretization for inerita
  real              :: urf           ! under-relaxation factor
!------------------------------------------------------------------------------!
!
!  Stress tensor on the face s:
!
!    T = mi * [    2*dU/dx     dU/dy+dV/dx   dU/dz+dW/dx  ]
!             [  dU/dy+dV/dx     2*dV/dy     dV/dz+dW/dy  ]
!             [  dU/dz+dW/dx   dV/dz+dW/dy     2*dW/dz    ]
!
!  The forces, acting on the cell face are:
!
!    Fx = T11*Sx + T12*Sy + T13*Sz
!    Fy = T21*Sx + T22*Sy + T23*Sz
!    Fz = T31*Sx + T32*Sy + T33*Sz
!
!  which could also be written in the compact form:
!
!    {F} = [T]{S}
!
!  U component:
!
!    Fx = Txx*Sx + Txy*Sy + Txz*Sz
!
!  V component:
!
!    Fy = Tyx*Sx + Tyy*Sy + Tyz*Sz
!
!  W component:
!
!    Fz = Tzx*Sx + Tzy*Sy + Tzz*Sz
!
!------------------------------------------------------------------------------!
!
!  The form of equations which I am solving:
!
!     /             /              /               /
!    |     du      |              |               |
!    | rho -- dV + | rho u u dS = | mu DIV u dS - | p dS
!    |     dt      |              |               |
!   /             /              /               /
!
!  Dimension of the system under consideration
!
!     [a]{u} = {b}   [kgm/s^2]   [N]
!
!  Dimensions of certain variables
!
!     a              [kg/s]
!     u, v, w        [m/s]
!     bu, bv, bw     [kgm/s^2]   [N]
!     p, pp          [kg/ms^2]   [N/m^2]
!     flux           [kg/s]
!     au*, av*, aw*  [kgm/s^2]   [N]
!     du*, dv*, dw*  [kgm/s^2]   [N]
!     cu*, cv*, cw*  [kgm/s^2]   [N]
!     Wall visc.      vis_wall [kg/(m*s)]
!==============================================================================!

  ! Take aliases
  a => sol % a
  b => sol % b

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0
  f_stress   = 0.0

  ! User function
  call User_Mod_Beginning_Of_Compute_Momentum(grid, dt, ini)

  ! Calculate velocity magnitude for normalization
  vel_max = 0.0
  do c = -grid % n_bnd_cells, grid % n_cells
    vel_max = max(vel_max, sqrt(u % n(c)**2 + v % n(c)**2 + w % n(c)**2))
  end do
  call Comm_Mod_Global_Max_Real(vel_max)

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      ui % oo(c) = ui % o(c)
      ui % o (c) = ui % n(c)
    end do
  end if

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Momentum(adv_scheme)
  call Control_Mod_Blending_Coefficient_For_Momentum(blend)

  ! Compute phimax and phimin
  if(adv_scheme .ne. CENTRAL) then
    call Calculate_Minimum_Maximum(grid, ui % n, ui_min, ui_max) ! or ui % o ?
    goto 1  ! why on Earth this?
  end if

  ! New values
1 do c = 1, grid % n_cells
    ui % a(c) = 0.0
    ui % c(c) = 0.0
  end do

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Central differencing
    uis = grid % f(s) * ui % n(c1) + (1.0 - grid % f(s)) * ui % n(c2)

    if(adv_scheme .ne. CENTRAL) then
      call Advection_Scheme(grid, uis, s, ui % n, ui_min, ui_max,  &
                            ui_i, ui_j, ui_k,                      &
                            di, dj, dk,                            &
                            adv_scheme, blend)
    end if

    ! Compute advection term
    if(c2  > 0) then
      ui % a(c1) = ui % a(c1) - flux(s) * uis
      ui % a(c2) = ui % a(c2) + flux(s) * uis
    else
      ui % a(c1) = ui % a(c1) - flux(s) * uis
    end if

    ! Store upwinded part of the advection term in "c"
    if(flux(s)  < 0) then   ! from c2 to c1
      ui % c(c1) = ui % c(c1) - flux(s)*ui % n(c2)
      if(c2  > 0) then
        ui % c(c2) = ui % c(c2) + flux(s)*ui % n(c2)
      end if
    else
      ui % c(c1) = ui % c(c1) - flux(s)*ui % n(c1)
      if(c2  > 0) then
        ui % c(c2) = ui % c(c2) + flux(s)*ui % n(c1)
      end if
    end if

  end do ! through faces

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + (ui % a(c) - ui % c(c))
  end do

  !---------------!
  !               !
  !   Diffusion   !
  !               !
  !---------------!

  ! Set c terms back to zero
  do c = 1, grid % n_cells
    ui % c(c) = 0.0
  end do

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    vis_eff = viscosity

    if(turbulence_model .ne. NONE .and.  &
       turbulence_model .ne. DNS) then
      vis_eff = vis_eff + fw(s)*vis_t(c1)+(1.0-fw(s))*vis_t(c2)
    end if

    if(turbulence_model .eq. HYBRID_LES_RANS) then
      vis_eff = fw(s)*vis_t_eff(c1)+(1.0-fw(s))*vis_t_eff(c2) + viscosity
    end if

    if(c2 < 0) then
      if( turbulence_model .eq. LES_SMAGORINSKY .or.  &
          turbulence_model .eq. LES_DYNAMIC     .or.  &
          turbulence_model .eq. LES_WALE) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          vis_eff = vis_wall(c1)
        end if
      end if
    end if

    if( turbulence_model .eq. K_EPS_ZETA_F     .or.  &
        turbulence_model .eq. HYBRID_LES_RANS  .or.  &
        turbulence_model .eq. K_EPS) then 
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          vis_eff = vis_wall(c1)
        end if
      end if
    end if

    ! Add influence of Re stresses
    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
      if(turbulence_model_variant .ne. STABILIZED) then
        if(ui % name .eq. 'U') then
          uu_f = fw(s) * uu % n(c1) + (1.0-fw(s)) * uu % n(c2)
          uv_f = fw(s) * uv % n(c1) + (1.0-fw(s)) * uv % n(c2)
          uw_f = fw(s) * uw % n(c1) + (1.0-fw(s)) * uw % n(c2)
          f_stress = - (  uu_f * grid % sx(s)  &
                        + uv_f * grid % sy(s)  &
                        + uw_f * grid % sz(s) )
        else if(ui % name .eq. 'V') then
          uv_f = fw(s) * uv % n(c1) + (1.0-fw(s)) * uv % n(c2)
          vv_f = fw(s) * vv % n(c1) + (1.0-fw(s)) * vv % n(c2)
          vw_f = fw(s) * vw % n(c1) + (1.0-fw(s)) * vw % n(c2)
          f_stress =  - (  uv_f * grid % sx(s)  &
                         + vv_f * grid % sy(s)  &
                         + vw_f * grid % sz(s) )
        else if(ui % name .eq. 'W') then
          uw_f = fw(s) * uw % n(c1) + (1.0-fw(s)) * uw % n(c2)
          vw_f = fw(s) * vw % n(c1) + (1.0-fw(s)) * vw % n(c2)
          ww_f = fw(s) * ww % n(c1) + (1.0-fw(s)) * ww % n(c2)
          f_stress =  - (  uw_f * grid % sx(s)  &
                         + vw_f * grid % sy(s)  &
                         + ww_f * grid % sz(s) )
        end if
      end if
    end if

    ui_i_f = fw(s)*ui_i(c1) + (1.0-fw(s))*ui_i(c2)
    ui_j_f = fw(s)*ui_j(c1) + (1.0-fw(s))*ui_j(c2)
    ui_k_f = fw(s)*ui_k(c1) + (1.0-fw(s))*ui_k(c2)
    uj_i_f = fw(s)*uj_i(c1) + (1.0-fw(s))*uj_i(c2)
    uk_i_f = fw(s)*uk_i(c1) + (1.0-fw(s))*uk_i(c2)

    ! Total (exact) viscous stress
    f_ex = vis_eff*(      2.0*ui_i_f  * si(s)      &
                    + (ui_j_f+uj_i_f) * sj(s)      &
                    + (ui_k_f+uk_i_f) * sk(s) )

    a0 = vis_eff * f_coef(s)

    ! Implicit viscous stress
    f_im = (   ui_i_f*di(s)                &
             + ui_j_f*dj(s)                &
             + ui_k_f*dk(s))*a0

    ! Cross diffusion part
    ui % c(c1) = ui % c(c1) + f_ex - f_im + f_stress * density
    if(c2  > 0) then
      ui % c(c2) = ui % c(c2) - f_ex + f_im - f_stress * density
    end if

    ! Compute the coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    a12 = a12  - min(flux(s), real(0.0))
    a21 = a21  + max(flux(s), real(0.0))

    ! Fill the system matrix
    if(c2 > 0) then
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % dia(c1))  = a % val(a % dia(c1))  + a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
      a % val(a % dia(c2))  = a % val(a % dia(c2))  + a21
    else if(c2  < 0) then
      ! Outflow is not included because it was causing problems
      if((Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)  .or.  &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL)    .or.  &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) .or.  &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL)) then
         ! (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1) = b(c1) + a12 * ui % n(c2)
      end if
    end if

  end do  ! through faces

  ! Here we clean up momentum from the false diffusion
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(turbulence_model_variant .ne. STABILIZED) then
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        vis_tS = (fw(s)*vis_t(c1)+(1.0-fw(s))*vis_t(c2))
        a0 = f_coef(s)*vis_tS
        vis_eff = vis_tS

        ui_i_f = fw(s) * ui_i(c1) + (1.0-fw(s)) * ui_i(c2)
        ui_j_f = fw(s) * ui_j(c1) + (1.0-fw(s)) * ui_j(c2)
        ui_k_f = fw(s) * ui_k(c1) + (1.0-fw(s)) * ui_k(c2)
        uj_i_f = fw(s) * uj_i(c1) + (1.0-fw(s)) * uj_i(c2)
        uk_i_f = fw(s) * uk_i(c1) + (1.0-fw(s)) * uk_i(c2)

        f_ex = vis_eff*( 2.0*ui_i_f         * si(s) &
                          + (ui_j_f+uj_i_f) * sj(s) &
                          + (ui_k_f+uk_i_f) * sk(s) )

        f_im = (  ui_i_f * di(s)  &
                + ui_j_f * dj(s)  &
                + ui_k_f * dk(s)) * vis_eff * f_coef(s)

        b(c1) = b(c1) - vis_eff * (ui % n(c2) -ui % n(c1)) * f_coef(s)  &
              - f_ex + f_im
        if(c2  > 0) then
          b(c2) = b(c2) + vis_eff * (ui % n(c2) -ui % n(c1)) * f_coef(s)  &
                + f_ex - f_im
        end if
      end do
    end if
  end if

  ! Fully implicit treatment for cross diffusion terms
  do c = 1, grid % n_cells
    b(c) = b(c) + ui % c(c)
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  call Control_Mod_Time_Integration_Scheme(td_scheme)

  ! Two time levels; linear interpolation
  if(td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = density * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c) = b(c) + a0 * ui % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = density * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0 * a0 * ui % o(c) - 0.5 * a0 * ui % oo(c)
    end do
  end if

  !---------------------------------!
  !                                 !
  !   Pressure term contributions   !
  !                                 !
  !---------------------------------!

  !--------------------------!
  !   Global pressure drop   !
  !--------------------------!
  if(ui % name .eq. 'U') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk % p_drop_x * grid % vol(c)
    end do
  else if(ui % name .eq. 'V') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk % p_drop_y * grid % vol(c)
    end do
  else if(ui % name .eq. 'W') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk % p_drop_z * grid % vol(c)
    end do
  end if

  !---------------------------------!
  !   Local pressure distribution   !
  !---------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) - h_i(c) * grid % vol(c)
  end do

  !--------------------!
  !   Buoyancy force   !
  !--------------------!
  if(buoyancy) then
    if(ui % name .eq. 'U') then
      do c = 1, grid % n_cells
        b(c) = b(c) - density * grav_x * (t % n(c) - t_ref)  &
             * grid % vol(c)
      end do
    else if(ui % name .eq. 'V') then
      do c = 1, grid % n_cells
        b(c) = b(c) - density * grav_y * (t % n(c) - t_ref)  &
             * grid % vol(c)
      end do
    else if(ui % name .eq. 'W') then
      do c = 1, grid % n_cells
        b(c) = b(c) - density * grav_z * (t % n(c) - t_ref)  &
             * grid % vol(c)
      end do
    end if
  end if

  !----------------------------------------!
  !   All other terms defined by the user  !
  !----------------------------------------!
  call User_Mod_Force(grid, ui, a, b)

  !-----------------------------------!
  !                                   !
  !   Solve the equations for u,v,w   !
  !                                   !
  !-----------------------------------!

  ! Set under-relaxation factor then overwrite with conrol file if specified
  urf = 0.8
  call Control_Mod_simple_Underrelaxation_For_Momentum(urf)

  do c = 1, grid % n_cells
    a % sav(c) = a % val(a % dia(c))
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - urf)*ui % n(c) / urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / urf
  end do

  ! Get solver tolerance
  call Control_Mod_Tolerance_For_Momentum_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set number of iterations then overwrite with conrol file if specified
  niter =  5
  call Control_Mod_Max_Iterations_For_Momentum_Solver(niter)

  call Bicg(sol,      &
            ui % n,   &
            b,        &
            precond,  &
            niter,    &
            tol,      &
            ini_res,  &
            ui % res, &
            norm = vel_max)

  if(ui % name .eq. 'U') then
    call Info_Mod_Iter_Fill_At(2, 1, ui % name, niter, ui % res)
  end if
  if(ui % name .eq. 'V') then
    call Info_Mod_Iter_Fill_At(2, 2, ui % name, niter, ui % res)
  end if
  if(ui % name .eq. 'W') then
    call Info_Mod_Iter_Fill_At(2, 3, ui % name, niter, ui % res)
  end if

  call Comm_Mod_Exchange_Real(grid, ui % n)

  ! User function
  call User_Mod_End_Of_Compute_Momentum(grid, dt, ini)

  end subroutine
