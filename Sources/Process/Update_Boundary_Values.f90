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
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!---------------------------------[Calling]------------------------------------!
  real :: Turbulent_Prandtl_Number
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
  real    :: qx, qy, qz, nx, ny, nz, s_tot, con_t, ebf, y_pl
  real    :: pr, beta, u_plus
!==============================================================================!

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
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.  &
          Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE ) then
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = w % n(c1)
        if(heat_transfer .eq. YES) t % n(c2) = t % n(c1)
      end if

      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY ) then
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = w % n(c1)
        if(heat_transfer .eq. YES) t % n(c2) = t % n(c1)
      end if

      ! Spalart Allmaras
      if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
         turbulence_model .eq. DES_SPALART) then
        if ( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE ) then
          vis % n(c2) = vis % n(c1)
        end if
      end if

      if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
         turbulence_model .eq. HANJALIC_JAKIRLIC) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          uu % n(c2) = 0.0
          vv % n(c2) = 0.0
          ww % n(c2) = 0.0
          uv % n(c2) = 0.0
          uw % n(c2) = 0.0
          vw % n(c2) = 0.0
          KIN % n(c2) = 0.0
          if(turbulence_model .eq. REYNOLDS_STRESS) f22% n(c2) = 0.0
        end if
      end if

      ! k-epsilon-zeta-f
      if(turbulence_model .eq. K_EPS_ZETA_F     .or.  &
         turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then
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

      if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
         turbulence_model .eq. HANJALIC_JAKIRLIC) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then
          kin % n(c2) = kin % n(c1)
          eps % n(c2) = eps % n(c1)
          uu % n(c2) = uu % n(c1)
          vv % n(c2) = vv % n(c1)
          ww % n(c2) = ww % n(c1)
          uv % n(c2) = uv % n(c1)
          uw % n(c2) = uw % n(c1)
          vw % n(c2) = vw % n(c1)
          if(turbulence_model .eq. REYNOLDS_STRESS) f22 % n(c2) = f22 % n(c1)
        end if
      end if

      ! Is this good in general case, when q <> 0 ??? Check it.
      if(heat_transfer .eq. YES) then
        call Control_Mod_Turbulent_Prandtl_Number(pr_t)
        if(turbulence_model .ne. LES .or.  &
           turbulence_model .ne. DNS) then
          pr_t = Turbulent_Prandtl_Number(grid, c1)
        end if
        s_tot = sqrt(  grid % sx(s)*grid % sx(s)  &
                     + grid % sy(s)*grid % sy(s)  &
                     + grid % sz(s)*grid % sz(s))
        nx = grid % sx(s)/s_tot
        ny = grid % sy(s)/s_tot
        nz = grid % sz(s)/s_tot
        qx = t % q(c2) * nx
        qy = t % q(c2) * ny
        qz = t % q(c2) * nz
        con_t = conductivity                 &
              + capacity*vis_t(c1) / pr_t
        if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
           turbulence_model .eq. K_EPS) then
          y_pl = max(c_mu25 * sqrt(kin % n(c1)) * grid % wall_dist(c1) / &
            viscosity, 0.12)
          u_plus = log(y_pl * e_log) / kappa + TINY
          pr = viscosity * capacity / conductivity
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)  &
                      * (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ebf = 0.01 * (pr*y_pl)**4.0          &
                     / (1.0 + 5.0 * pr**3 * y_pl) + TINY
          con_wall(c1) = y_pl * viscosity * capacity   &
                       / (y_pl * pr * exp(-1.0 * ebf)  &
                       + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + pr_t / capacity  &
                      * (  qx * grid % dx(s)         &
                         + qy * grid % dy(s)         &
                         + qz * grid % dz(s))        &
                      / con_wall(c1)
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * con_t  &
                      / (  nx * grid % dx(s)               &
                         + ny * grid % dy(s)               &
                         + nz * grid % dz(s) )
          end if
        else
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + pr_t / (con_t + TINY)   &
                      * (  qx * grid % dx(s)                &
                         + qy * grid % dy(s)                &
                         + qz * grid % dz(s) )
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * con_t  &
                      / (  nx * grid % dx(s)               &
                         + ny * grid % dy(s)               &
                         + nz * grid % dz(s) )
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

        if(heat_transfer .eq. YES)  &
          t % n(c2) = t % n(grid % bnd_cond % copy_c(c2))

        if(turbulence_model .eq. SPALART_ALLMARAS .or.   &
           turbulence_model .eq. DES_SPALART) &
          vis % n(c2) = vis % n(grid % bnd_cond % copy_c(c2))

        if(turbulence_model .eq. K_EPS_ZETA_F     .or.  &
           turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
          kin  % n(c2) = kin  % n(grid % bnd_cond % copy_c(c2))
          eps  % n(c2) = eps  % n(grid % bnd_cond % copy_c(c2))
          zeta % n(c2) = zeta % n(grid % bnd_cond % copy_c(c2))
          f22  % n(c2) = f22  % n(grid % bnd_cond % copy_c(c2))
        end if

        if(turbulence_model .eq. K_EPS) then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
        end if

        if(turbulence_model .eq. REYNOLDS_STRESS .or.  &
           turbulence_model .eq. HANJALIC_JAKIRLIC) then
          kin % n(c2) = kin % n(grid % bnd_cond % copy_c(c2))
          eps % n(c2) = eps % n(grid % bnd_cond % copy_c(c2))
          uu % n(c2)  = uu  % n(grid % bnd_cond % copy_c(c2))
          vv % n(c2)  = vv  % n(grid % bnd_cond % copy_c(c2))
          ww % n(c2)  = ww  % n(grid % bnd_cond % copy_c(c2))
          uv % n(c2)  = uv  % n(grid % bnd_cond % copy_c(c2))
          uw % n(c2)  = uw  % n(grid % bnd_cond % copy_c(c2))
          vw % n(c2)  = vw  % n(grid % bnd_cond % copy_c(c2))
          f22 % n(c2) = f22 % n(grid % bnd_cond % copy_c(c2))
        end if
      end if
    end if
  end do

  end subroutine
