!==============================================================================!
  subroutine Compute_Momentum(grid, dt, ini, ui,  &
                              ui_i, ui_j, ui_k,   &
                              si, sj, sk,         &
                              di, dj, dk,         &
                              Hi, uj_i, uk_i)
!------------------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Comm_Mod
  use Var_Mod
  use Grid_Mod
  use Bulk_Mod
  use Info_Mod
  use Numerics_Mod
  use Solvers_Mod, only: Bicg, Cg, Cgs
  use Control_Mod
  use User_Mod
  use Work_Mod,    only: ui_min => r_cell_01,  &
                         ui_max => r_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
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
  real            :: Hi  (-grid % n_bnd_cells:grid % n_cells),  &
                     uj_i(-grid % n_bnd_cells:grid % n_cells),  &
                     uk_i(-grid % n_bnd_cells:grid % n_cells)
  real            :: uu_f, vv_f, ww_f, uv_f, uw_f, vw_f
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter, mat
  real              :: f_ex, f_im, f_stress
  real              :: uis
  real              :: a0, a12, a21
  real              :: ini_res, tol
  real              :: vis_eff, vis_tS
  real              :: ui_i_f,ui_j_f,ui_k_f,uj_i_f,uk_i_f
  character(len=80) :: coupling
  character(len=80) :: precond
  integer           :: adv_scheme    ! space disretization of advection (scheme)
  real              :: blend         ! blending coeff (1.0 central; 0.0 upwind)
  integer           :: td_inertia    ! time-disretization for inerita
  integer           :: td_advection  ! time-disretization for advection
  integer           :: td_diffusion  ! time-disretization for diffusion
  integer           :: td_cross_diff ! time-disretization for cross-difusion
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

  ! User function
  call User_Mod_Beginning_Of_Compute_Momentum(grid, dt, ini)

  b = 0.0
  a % val = 0.0
  f_stress = 0.0

  ! This is important for "copy" boundary conditions. Find out why !
  do c = -grid % n_bnd_cells, -1
    a % bou(c) = 0.0
  end do

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
      ui % oo(c)   = ui % o(c)
      ui % o (c)   = ui % n(c)
      ui % a_oo(c) = ui % a_o(c)
      ui % a_o (c) = 0.0
      ui % d_oo(c) = ui % d_o(c)
      ui % d_o (c) = 0.0
      ui % c_oo(c) = ui % c_o(c)
      ui % c_o (c) = ui % c(c)
    end do
  end if

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Momentum(adv_scheme)
  call Control_Mod_Blending_Coefficient_Momentum(blend)

  ! Compute phimax and phimin
  do mat = 1, grid % n_materials
    if(adv_scheme .ne. CENTRAL) then
      call Calculate_Minimum_Maximum(grid, ui % n, ui_min, ui_max) ! or ui % o ?
      goto 1  ! why on Earth this?
    end if
  end do

  ! New values
1 do c = 1, grid % n_cells
    ui % a(c)    = 0.0
    ui % c(c)    = 0.0
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
    if(ini .eq. 1) then
      if(c2  > 0) then
        ui % a_o(c1) = ui % a_o(c1) - flux(s) * uis
        ui % a_o(c2) = ui % a_o(c2) + flux(s) * uis
      else
        ui % a_o(c1) = ui % a_o(c1) - flux(s) * uis
      end if
    end if

    if(c2  > 0) then
      ui % a(c1) = ui % a(c1) - flux(s) * uis
      ui % a(c2) = ui % a(c2) + flux(s) * uis
    else
      ui % a(c1) = ui % a(c1) - flux(s) * uis
    end if

    ! Store upwinded part of the advection term in "c"
    if(coupling .ne. 'PROJECTION') then
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
    end if
  end do ! through faces

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  ! Adams-Bashforth scheeme for convective fluxes
  if(td_advection .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c) = b(c) + (1.5*ui % a_o(c) - 0.5*ui % a_oo(c) - ui % c(c))
    end do
  end if

  ! Crank-Nicholson scheeme for convective fluxes
  if(td_advection .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c) = b(c) + (0.5 * ( ui % a(c) + ui % a_o(c) ) - ui % c(c))
    end do
  end if

  ! Fully implicit treatment of convective fluxes
  if(td_advection .eq. FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      b(c) = b(c) + (ui % a(c) - ui % c(c))
    end do
  end if

  ! New values
  do c = 1, grid % n_cells
    ui % c(c) = 0.0
  end do

  !---------------!
  !               !
  !   Diffusion   !
  !               !
  !---------------!

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

    if(turbulence_model .eq. K_EPS_ZETA_F .and.  &
       turbulence_statistics .eq. YES) then
      vis_eff = fw(s)*vis_t_eff(c1)+(1.0-fw(s))*vis_t_eff(c2) + viscosity
    end if

    if(c2 < 0) then
      if( turbulence_model .eq. SMAGORINSKY .or.  &
          turbulence_model .eq. DYNAMIC     .or.  &
          turbulence_model .eq. WALE) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          vis_eff = vis_wall(c1)
        end if
      end if
    end if

    if( turbulence_model .eq. K_EPS_ZETA_F     .or.  &
       (turbulence_model .eq. K_EPS .and.            &
        turbulence_wall_treatment .eq. HIGH_RE) ) then
      if(c2 < 0) then
        if (Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            vis_eff = vis_wall(c1)
          end if
        end if
      end if
    end if

    ! Add influence of Re stresses
    if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
       turbulence_model .eq. HANJALIC_JAKIRLIC) then
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

    ! total (exact) viscous stress
    f_ex = vis_eff*(      2.0*ui_i_f  * si(s)      &
                    + (ui_j_f+uj_i_f) * sj(s)      &
                    + (ui_k_f+uk_i_f) * sk(s) )

    a0 = vis_eff * f_coef(s)

    ! Implicit viscous stress
    ! this is a very crude approximation: f_coef is not
    ! corrected at interface between materials
    f_im = (   ui_i_f*di(s)                &
             + ui_j_f*dj(s)                &
             + ui_k_f*dk(s))*a0

    ! Straight diffusion part
    if(ini .eq. 1) then
      if(c2  > 0) then
        ui % d_o(c1) = ui % d_o(c1) + (ui % n(c2)-ui % n(c1))*a0
        ui % d_o(c2) = ui % d_o(c2) - (ui % n(c2)-ui % n(c1))*a0
      else
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. SYMMETRY) then
          ui % d_o(c1) = ui % d_o(c1) + (ui % n(c2)-ui % n(c1))*a0
        end if
      end if
    end if

    ! Cross diffusion part
    ui % c(c1) = ui % c(c1) + f_ex - f_im + f_stress
    if(c2  > 0) then
      ui % c(c2) = ui % c(c2) - f_ex + f_im - f_stress
    end if

    ! Compute the coefficients for the sysytem matrix
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

      if(coupling .ne. 'PROJECTION') then
        a12 = a12  - min(flux(s), real(0.0))
        a21 = a21  + max(flux(s), real(0.0))
      end if

      ! Fill the system matrix
      if(c2  > 0) then
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
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then
          a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
          a % bou(c2) = -a12  ! cool parallel stuff
        end if
      end if
    end if
  end do  ! through faces

  !-----------------------------!
  !   Temporal discretization   !
  !-----------------------------!

  !---------------------------------!
  !
  ! Add Re stress influence on momentum
  ! This is an alternative way to implement RSM
  !
  !  if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
  !     turbulence_model .eq. HANJALIC_JAKIRLIC) then
  !    if(ui % name .eq. 'U') then
  !      call GraPhi(uu%n,1,VAR2x,.TRUE.)
  !      call GraPhi(uv%n,2,VAR2y,.TRUE.)
  !      call GraPhi(uw%n,3,VAR2z,.TRUE.)
  !      do c = 1, grid % n_cells
  !        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*grid % vol(c)
  !      end do
  !    else if(ui % name .eq. 'V') then
  !      call GraPhi(uv%n,1,VAR2x,.TRUE.)
  !      call GraPhi(vv%n,2,VAR2y,.TRUE.)
  !      call GraPhi(vw%n,3,VAR2z,.TRUE.)
  !      do c = 1, grid % n_cells
  !        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*grid % vol(c)
  !      end do
  !    else if(ui % name .eq. 'W') then
  !      call GraPhi(uw%n,1,VAR2x,.TRUE.)
  !      call GraPhi(vw%n,2,VAR2y,.TRUE.)
  !      call GraPhi(ww%n,3,VAR2z,.TRUE.)
  !      do c = 1, grid % n_cells
  !        b(c) = b(c) - (VAR2x(c)+VAR2y(c)+VAR2z(c))*grid % vol(c)
  !      end do
  !    end if

  ! Here we clean up momentum from the false diffusion
  if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
     turbulence_model .eq. HANJALIC_JAKIRLIC) then
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

  ! Adams-Bashfort scheeme for diffusion fluxes
  if(td_diffusion .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * ui % d_o(c) - 0.5 * ui % d_oo(c)
    end do
  end if

  ! Crank-Nicholson scheme for difusive terms
  if(td_diffusion .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * ui % d_o(c)
    end do
  end if

  ! Fully implicit treatment for difusive terms
  ! is handled via the linear system of equations

  ! Adams-Bashfort scheeme for cross diffusion
  if(td_cross_diff .eq. ADAMS_BASHFORTH) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 1.5 * ui % c_o(c) - 0.5 * ui % c_oo(c)
    end do
  end if

  ! Crank-Nicholson scheme for cross difusive terms
  if(td_cross_diff .eq. CRANK_NICOLSON) then
    do c = 1, grid % n_cells
      b(c) = b(c) + 0.5 * ui % c(c) + 0.5 * ui % c_o(c)
    end do
  end if

  ! Fully implicit treatment for cross difusive terms
  if(td_cross_diff .eq. FULLY_IMPLICIT) then
    do c = 1, grid % n_cells
      b(c) = b(c) + ui % c(c)
    end do
  end if

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  ! Two time levels; linear interpolation
  if(td_inertia .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = density * grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c) = b(c) + a0 * ui % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_inertia .eq. PARABOLIC) then
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
      b(c) = b(c) + bulk(grid % material(c)) % p_drop_x * grid % vol(c)
    end do
  else if(ui % name .eq. 'V') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk(grid % material(c)) % p_drop_y * grid % vol(c)
    end do
  else if(ui % name .eq. 'W') then
    do c = 1, grid % n_cells
      b(c) = b(c) + bulk(grid % material(c)) % p_drop_z * grid % vol(c)
    end do
  end if

  !---------------------------------!
  !   Local pressure distribution   !
  !---------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) - Hi(c) * grid % vol(c)
  end do

  !----------------------------------------!
  !   All other terms defined by the user  !
  !----------------------------------------!
  call User_Mod_Force(grid, ui, a, b)

  !-----------------------------------!
  !                                   !
  !   Solve the equations for u,v,w   !
  !                                   !
  !-----------------------------------!

  ! Type of coupling is important
  call Control_Mod_Pressure_Momentum_Coupling(coupling)

  ! Set under-relaxation factor
  urf = 1.0
  if(coupling .eq. 'SIMPLE')  &
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

  ! Set number of solver iterations on coupling method
  if(coupling .eq. 'PROJECTION') niter = 10
  if(coupling .eq. 'SIMPLE')     niter =  5

  ! Over-ride if specified in control file
  call Control_Mod_Max_Iterations_For_Momentum_Solver(niter)

  call Cg(a, ui % n, b, precond, niter, tol, ini_res, ui % res)

  if(ui % name .eq. 'U') then
    call Info_Mod_Iter_Fill_At(2, 1, ui % name, niter, ui % res)
  end if
  if(ui % name .eq. 'V') then
    call Info_Mod_Iter_Fill_At(2, 2, ui % name, niter, ui % res)
  end if
  if(ui % name .eq. 'W') then
    call Info_Mod_Iter_Fill_At(2, 3, ui % name, niter, ui % res)
  end if

  call Comm_Mod_Exchange(grid, ui % n)

  ! User function
  call User_Mod_End_Of_Compute_Momentum(grid, dt, ini)

  end subroutine
