!==============================================================================!
  subroutine Update_Boundary_Values(flow, turb, mult)
!------------------------------------------------------------------------------!
!   Update variables on the boundaries (boundary cells) where needed.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Field_Mod,   only: Field_Type, heat_transfer
  use Turb_Mod
  use Grid_Mod
  use Control_Mod
  use Multiphase_Mod, only: Multiphase_Type, Multiphase_Mod_Alias_Vof,  &
                            multiphase_model, VOLUME_OF_FLUID
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, phi, vof
  type(Var_Type),  pointer :: kin, eps, zeta, f22, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  integer                  :: c1, c2, s, sc
  real                     :: qx, qy, qz, nx, ny, nz
  real                     :: kin_vis, u_tau
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  vis  => turb % vis
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Field_Mod_Alias_Energy     (flow, t)
  call Multiphase_Mod_Alias_Vof   (mult, vof)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_T2          (turb, t2)

  if(heat_transfer) then
    flow % heat        = 0.0     ! [W]
    flow % heat_flux   = 0.0     ! [W/m^2]
    flow % heated_area = 0.0
    ! Take (default) turbulent Prandtl number from control file
    call Control_Mod_Turbulent_Prandtl_Number(pr_t)
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    !------------------------------------------------!
    !   Outflow (and inflow, if needed) boundaries   !
    !------------------------------------------------!

    ! On the boundary perform the extrapolation
    if(c2 < 0) then

      kin_vis = flow % viscosity(c1) / flow % density(c1)
      ! Extrapolate velocities on the outflow boundary
      ! SYMMETRY is intentionally not treated here because I wanted to
      ! be sure that is handled only via graPHI and NewUVW functions)
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.     &
          Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.    &
          Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY ) then
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = w % n(c1)
        if(heat_transfer) t % n(c2) = t % n(c1)
        if(multiphase_model .eq. VOLUME_OF_FLUID) vof % n(c2) = vof % n(c1)
      end if

      ! Spalart Allmaras
      if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
         turbulence_model .eq. DES_SPALART) then
        if ( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
          vis % n(c2) = vis % n(c1)
        end if
      end if

      ! Reynolds stress models
      if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
         turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          uu  % n(c2) = 0.0
          vv  % n(c2) = 0.0
          ww  % n(c2) = 0.0
          uv  % n(c2) = 0.0
          uw  % n(c2) = 0.0
          vw  % n(c2) = 0.0
          kin % n(c2) = 0.0
          u_tau = kin_vis * sqrt(u % n(c1)**2 + v % n(c1)**2 + w % n(c1)**2)  &
                / grid % wall_dist(c1)
          turb % y_plus(c1) = Y_Plus_Low_Re(turb, u_tau,           &
                                            grid % wall_dist(c1),  &
                                            kin_vis)
          if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) f22 % n(c2) = 0.0
        end if
      end if

      ! k-epsilon-zeta-f
      if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
         turbulence_model .eq. HYBRID_LES_RANS) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.   &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.   &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.   &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
          kin  % n(c2) = kin  % n(c1)
          eps  % n(c2) = eps  % n(c1)
          zeta % n(c2) = zeta % n(c1)
          f22  % n(c2) = f22  % n(c1)
          if(heat_transfer) then
            t2  % n(c2) = t2  % n(c1)
          end if
        end if
      end if

      ! k-epsilon
      if(turbulence_model .eq. K_EPS) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
          if(heat_transfer) then
            t2  % n(c2) = t2  % n(c1)
          end if 
        end if
      end if

      if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
         turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then
          uu  % n(c2) = uu  % n(c1)
          vv  % n(c2) = vv  % n(c1)
          ww  % n(c2) = ww  % n(c1)
          uv  % n(c2) = uv  % n(c1)
          uw  % n(c2) = uw  % n(c1)
          vw  % n(c2) = vw  % n(c1)
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
          if(turbulence_model .eq. RSM_MANCEAU_HANJALIC)  &
            f22 % n(c2) = f22 % n(c1)
        end if
      end if

      ! Is this good in general case, when q <> 0 ??? Check it.
      if(heat_transfer) then

        ! If not DNS or LES, compute Prandtl number 
        if(turbulence_model .ne. LES_SMAGORINSKY     .and.  &
           turbulence_model .ne. LES_DYNAMIC         .and.  &
           turbulence_model .ne. LES_WALE            .and.  &
           turbulence_model .ne. HYBRID_LES_PRANDTL  .and.  &
           turbulence_model .ne. NONE                .and.  &
           turbulence_model .ne. DNS) then
          pr_t = Turb_Mod_Prandtl_Number(turb, c1)
        end if

        nx = grid % sx(s) / grid % s(s)
        ny = grid % sy(s) / grid % s(s)
        nz = grid % sz(s) / grid % s(s)
        qx = t % q(c2) * nx
        qy = t % q(c2) * ny
        qz = t % q(c2) * nz

        ! Wall temperature or heat fluxes for k-eps-zeta-f
        ! and high-re k-eps models. 
        if(turbulence_model .eq. K_EPS_ZETA_F    .or.  &
           turbulence_model .eq. HYBRID_LES_RANS .or.  &
           turbulence_model .eq. K_EPS) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * grid % wall_dist(c1)  &
                      / (turb % con_w(c1) + TINY)
            flow % heat = flow % heat + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) then
              flow % heated_area = flow % heated_area + grid % s(s)
            end if
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * turb % con_w(c1)  &
                      / grid % wall_dist(c1)
            flow % heat = flow % heat + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) then
              flow % heated_area = flow % heated_area + grid % s(s)
            end if
          end if

        ! Wall temperature or heat fluxes for other trubulence models
        else
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * grid % wall_dist(c1)  &
                      / flow % conductivity(c1)
            flow % heat = flow % heat + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) then
              flow % heated_area = flow % heated_area + grid % s(s)
            end if
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * flow % conductivity(c1)  &
                      / grid % wall_dist(c1)
            flow % heat = flow % heat + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) then
              flow % heated_area = flow % heated_area + grid % s(s)
            end if
          end if
        end if
      end if ! heat_transfer

    end if ! c2 < 0
  end do

  ! Integrate (summ) heated area, and heat up
  if(heat_transfer) then
    call Comm_Mod_Global_Sum_Real(flow % heat)
    call Comm_Mod_Global_Sum_Real(flow % heated_area)
    flow % heat_flux = flow % heat / max(flow % heated_area, TINY)
  end if

  !-------------!
  !   Scalars   !
  !-------------!
  do sc = 1, flow % n_scalars
    phi => flow % scalar(sc)

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if (c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          phi % n(c2) = phi % n(c1) + phi % q(c2) * grid % wall_dist(c1)  &
                      / flow % diffusivity
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
          phi % q(c2) = (phi % n(c2) - phi % n(c1)) * flow % diffusivity &
                      / grid % wall_dist(c1)
        end if ! WALL or WALLFL

      end if ! c2 < 0
    end do ! s = 1, grid % n_faces
  end do ! sc = 1, flow % n_scalars

  end subroutine
