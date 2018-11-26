!==============================================================================!
  subroutine Update_Boundary_Values(grid)
!------------------------------------------------------------------------------!
!   Update variables on the boundaries (boundary cells) where needed.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Grid_Mod
  use Comm_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!---------------------------------[Calling]------------------------------------!
  real :: Turbulent_Prandtl_Number
  real :: Y_Plus_Low_Re
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
  real    :: qx, qy, qz, nx, ny, nz, con_t, ebf 
  real    :: pr, beta, u_plus, heated_area, y_pl, kin_vis
!==============================================================================!

  heat_flux   = 0.0
  heated_area = 0.0
  kin_vis     = viscosity / density

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    !------------------------------------------------!
    !   Outflow (and inflow, if needed) boundaries   !
    !------------------------------------------------!

    ! On the boundary perform the extrapolation
    if(c2  < 0) then

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
          u_tau(c1)   = kin_vis * sqrt(u % n(c1)**2      &
                        + v % n(c1)**2 + w % n(c1)**2)   &
                        / grid % wall_dist(c1)
          y_plus(c1)  = Y_Plus_Low_Re(u_tau(c1),         &
                        grid % wall_dist(c1), kin_vis) 
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
        end if

        !  if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        !    f22 % n(c2) = f22 % n(c1)
        !  end if
      end if

      ! k-epsilon
      if(turbulence_model .eq. K_EPS) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
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
          if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) f22 % n(c2) = f22 % n(c1)
        end if
      end if

      ! Is this good in general case, when q <> 0 ??? Check it.
      if(heat_transfer) then

        ! Take (default) turbulent Prandtl number from control file
        call Control_Mod_Turbulent_Prandtl_Number(pr_t)

        ! If not DNS or LES, compute Prandtl number 
        if(turbulence_model .ne. LES_SMAGORINSKY .and.  &
           turbulence_model .ne. LES_DYNAMIC     .and.  &
           turbulence_model .ne. LES_WALE        .and.  &
           turbulence_model .ne. NONE            .and.  &
           turbulence_model .ne. DNS) then
          pr_t = Turbulent_Prandtl_Number(grid, c1)
        end if

        nx = grid % sx(s) / grid % s(s)
        ny = grid % sy(s) / grid % s(s)
        nz = grid % sz(s) / grid % s(s)
        qx = t % q(c2) * nx
        qy = t % q(c2) * ny
        qz = t % q(c2) * nz

        ! Turbulent conductivity from Reynolds analogy
        if(turbulence_model .ne. NONE .and.    &
           turbulence_model .ne. DNS) then
          con_t = conductivity                 &
                + capacity*vis_t(c1) / pr_t
        else
          con_t = conductivity
        end if

        ! Wall temperature or heat fluxes for k-eps-zeta-f
        ! and high-re k-eps models
        if(turbulence_model .eq. K_EPS_ZETA_F    .or.  &
           turbulence_model .eq. HYBRID_LES_RANS .or.  &
           turbulence_model .eq. K_EPS) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * grid % wall_dist(c1)  &
                      / (con_wall(c1) + TINY)
            heat_flux = heat_flux + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) heated_area = heated_area + grid % s(s)
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * con_wall(c1)  &
                      / grid % wall_dist(c1)
            heat_flux = heat_flux + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) heated_area = heated_area + grid % s(s)
          end if

        ! Wall temperature or heat fluxes for other trubulence models
        else
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * grid % wall_dist(c1)  &
                       /conductivity
            heat_flux = heat_flux + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) heated_area = heated_area + grid % s(s)
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * conductivity  &
                      / grid % wall_dist(c1)
            heat_flux = heat_flux + t % q(c2) * grid % s(s)
            if(abs(t % q(c2)) > TINY) heated_area = heated_area + grid % s(s)
          end if
        end if
      end if

      !---------------------!
      !   Copy boundaries   !
      !---------------------!
      if(grid % bnd_cond % copy_c(c2) .ne. 0) then
        u % n(c2) = u % n(grid % bnd_cond % copy_c(c2))
        v % n(c2) = v % n(grid % bnd_cond % copy_c(c2))
        w % n(c2) = w % n(grid % bnd_cond % copy_c(c2))

        if(heat_transfer)  &
          t % n(c2) = t % n(grid % bnd_cond % copy_c(c2))

        if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
           turbulence_model .eq. DES_SPALART)           &
          vis % n(c2) = vis % n(grid % bnd_cond % copy_c(c2))

        if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
           turbulence_model .eq. HYBRID_LES_RANS) then
          kin  % n(c2) = kin  % n(grid % bnd_cond % copy_c(c2))
          eps  % n(c2) = eps  % n(grid % bnd_cond % copy_c(c2))
          zeta % n(c2) = zeta % n(grid % bnd_cond % copy_c(c2))
          f22  % n(c2) = f22  % n(grid % bnd_cond % copy_c(c2))
        end if

        if(turbulence_model .eq. K_EPS) then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
        end if

        if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
           turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
          uu  % n(c2) = uu  % n(grid % bnd_cond % copy_c(c2))
          vv  % n(c2) = vv  % n(grid % bnd_cond % copy_c(c2))
          ww  % n(c2) = ww  % n(grid % bnd_cond % copy_c(c2))
          uv  % n(c2) = uv  % n(grid % bnd_cond % copy_c(c2))
          uw  % n(c2) = uw  % n(grid % bnd_cond % copy_c(c2))
          vw  % n(c2) = vw  % n(grid % bnd_cond % copy_c(c2))
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
          if(turbulence_model .eq. RSM_MANCEAU_HANJALIC)  &
            f22 % n(c2) = f22 % n(grid % bnd_cond % copy_c(c2))
        end if
      end if
    end if
  end do

  if(heat_transfer) then
    call Comm_Mod_Global_Sum_Real(heat_flux)
    call Comm_Mod_Global_Sum_Real(heated_area)
    heat_flux = heat_flux / (heated_area + TINY)
    heat      = heat_flux * heated_area
  end if

  end subroutine
