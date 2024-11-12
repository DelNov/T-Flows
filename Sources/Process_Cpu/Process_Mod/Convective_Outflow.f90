!==============================================================================!
  subroutine Convective_Outflow(Process, Flow, Turb, Vof)
!------------------------------------------------------------------------------!
!>  This subroutine implements convective outlflow condition from A. Bottaro.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Sets up pointers and aliases for grid, flow field,       !
!     turbulence, and multiphase variables. Prepares for boundary value        !
!     extrapolation based on the simulation timestep.                          !
!   * Momentum extrapolation: Depending on the simulation timestep, updates    !
!     the momentum variables (u, v, w) at outflow boundaries using either a    !
!     direct extrapolation or a convective method.                             !
!   * Turbulence variables: Updates turbulence quantities like kinetic energy, !
!     dissipation, and stress tensors at convective outflow boundaries based   !
!     on the turbulence model (K-epsilon, Spalart-Allmaras, RSM, etc.).        !
!   * Scalar variables: For simulations with scalar transport, extrapolates    !
!     scalar variables at outflow boundaries using a similar approach as       !
!     momentum variables.                                                      !
!   * Energy variables: In simulations involving heat transfer, updates        !
!     temperature at outflow boundaries, considering the flow's bulk velocity  !
!     and temperature gradients.                                               !
!   * Boundary condition handling: Applies different extrapolation strategies  !
!     depending on whether the simulation timestep is before or after a        !
!     specified threshold (BEGIN). This allows for finer control over boundary !
!     value updates during the simulation.                                     !
!   * Performance monitoring: Monitors the subroutine's performance, aiding    !
!     in optimizing simulation efficiency and accuracy.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Turb_Type),     target :: Turb     !! turbulence object
  type(Vof_Type),      target :: Vof      !! VOF object
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: BEGIN = 12
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t, phi, vis
  type(Var_Type),  pointer :: kin, eps, zeta, f22, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Face_Type), pointer :: v_flux
  integer                  :: c1, c2, s, sc, reg
  real                     :: nx, ny, nz, bulk_vel, phi_n, dt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  v_flux => Flow % v_flux
  dt     =  Flow % dt
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_T2          (t2)
  call Turb % Alias_Vis         (vis)

  !------------------------------------------------!
  !   Compute bulk velocity via a user function.   !
  !- - - - - - - - - - - - - - - - - - - - - - - - !
  !   Although the default version of the called   !
  !   function User_Mod_Bulk_Velocity should do    !
  !   a decent job in most cases, there could be   !
  !   some in which it is trickier to define.      !
  !   Cases with multiple outflows come to mind,   !
  !   or internal sinks or sources,                !
  !------------------------------------------------!
  call User_Mod_Bulk_Velocity(Flow, bulk_vel)

  !------------------------!
  !                        !
  !   Momentum equations   !
  !                        !
  !------------------------!

  if(Time % Curr_Dt() > BEGIN) then

    call Flow % Grad_Variable(Flow % u)
    call Flow % Grad_Variable(Flow % v)
    call Flow % Grad_Variable(Flow % w)

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. CONVECT) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          call Grid % Face_Normal(s, nx, ny, nz)

          phi_n = u % x(c1) * nx + u % y(c1) * ny + u % z(c1) * nz
          u % n(c2) = u % n(c2) - bulk_vel * phi_n * dt

          phi_n = v % x(c1) * nx + v % y(c1) * ny + v % z(c1) * nz
          v % n(c2) = v % n(c2) - bulk_vel * phi_n * dt

          phi_n = w % x(c1) * nx + w % y(c1) * ny + w % z(c1) * nz
          w % n(c2) = w % n(c2) - bulk_vel * phi_n * dt
        end do    ! face
      end if      ! boundary condition
    end do        ! region

  else      ! Time % Curr_Dt() <= BEGIN

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. CONVECT) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          u % n(c2) = u % n(c1)
          v % n(c2) = v % n(c1)
          w % n(c2) = w % n(c1)
        end do    ! face
      end if      ! boundary condition
    end do        ! region

  end if    ! Time % Curr_Dt() > BEGIN

  !--------------------------!
  !                          !
  !   Turbulent quantities   !
  !                          !
  !--------------------------!

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(Turb % model .eq. K_EPS) then

    if(Time % Curr_Dt() > BEGIN) then

      call Flow % Grad_Variable(kin)
      call Flow % Grad_Variable(eps)
      if(Flow % heat_transfer) then
        call Flow % Grad_Variable(t2)
      end if

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            call Grid % Face_Normal(s, nx, ny, nz)

            phi_n = kin % x(c1) * nx + kin % y(c1) * ny + kin % z(c1) * nz
            kin % n(c2) = kin % n(c2) - bulk_vel * phi_n * dt

            phi_n = eps % x(c1) * nx + eps % y(c1) * ny + eps % z(c1) * nz
            eps % n(c2) = eps % n(c2) - bulk_vel * phi_n * dt

            if(Flow % heat_transfer) then
              phi_n = t2 % x(c1) * nx + t2 % y(c1) * ny + t2 % z(c1) * nz
              t2 % n(c2) = t2 % n(c2) - bulk_vel * phi_n * dt
            end if
          end do    ! face
        end if      ! boundary condition
      end do        ! region

    else      ! Time % Curr_Dt() <= BEGIN

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            kin % n(c2) = kin % n(c1)
            eps % n(c2) = eps % n(c1)
            if(Flow % heat_transfer) then
              t2 % n(c2) = t2 % n(c1)
            end if
          end do    ! face
        end if      ! boundary condition
      end do        ! region

    end if    ! Time % Curr_Dt() > BEGIN

  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then

    if(Time % Curr_Dt() > BEGIN) then

      call Flow % Grad_Variable(kin)
      call Flow % Grad_Variable(eps)
      call Flow % Grad_Variable(f22)
      call Flow % Grad_Variable(zeta)
      if(Flow % heat_transfer) then
        call Flow % Grad_Variable(t2)
      end if

      ! On the boundary perform the extrapolation
      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            call Grid % Face_Normal(s, nx, ny, nz)

            phi_n = kin % x(c1) * nx + kin % y(c1) * ny + kin % z(c1) * nz
            kin % n(c2) = kin % n(c2) - bulk_vel * phi_n * dt

            phi_n = eps % x(c1) * nx + eps % y(c1) * ny + eps % z(c1) * nz
            eps % n(c2) = eps % n(c2) - bulk_vel * phi_n * dt

            phi_n = f22 % x(c1) * nx + f22 % y(c1) * ny + f22 % z(c1) * nz
            f22 % n(c2) = f22 % n(c2) - bulk_vel * phi_n * dt

            phi_n = zeta % x(c1) * nx + zeta % y(c1) * ny + zeta % z(c1) * nz
            zeta % n(c2) = zeta % n(c2) - bulk_vel * phi_n * dt

            if(Flow % heat_transfer) then
              phi_n = t2 % x(c1) * nx + t2 % y(c1) * ny + t2 % z(c1) * nz
              t2 % n(c2) = t2 % n(c2) - bulk_vel * phi_n * dt
            end if
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    else      ! Time % Curr_Dt() <= BEGIN

      ! On the boundary perform the extrapolation
      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            kin % n(c2)  = kin % n(c1)
            eps % n(c2)  = eps % n(c1)
            f22 % n(c2)  = f22 % n(c1)
            zeta % n(c2) = zeta % n(c1)
            if(Flow % heat_transfer) then
              t2 % n(c2) = t2 % n(c1)
            end if
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    end if    ! Time % Curr_Dt() > BEGIN

  end if

  !----------------------!
  !   Spalart-Allmaras   !
  !----------------------!
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then

    if(Time % Curr_Dt() > BEGIN) then

      call Flow % Grad_Variable(vis)

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            call Grid % Face_Normal(s, nx, ny, nz)

            phi_n = vis % x(c1) * nx + vis % y(c1) * ny + vis % z(c1) * nz
            vis % n(c2) = vis % n(c2) - bulk_vel * phi_n * dt
          end do    ! face
        end if      ! boundary condition
      end do        ! region

    else      ! Time % Curr_Dt() <= BEGIN

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            vis % n(c2) = vis % n(c1)
          end do    ! face
        end if      ! boundary condition
      end do        ! region

    end if    ! Time % Curr_Dt() > BEGIN

  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!

  if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    if(Time % Curr_Dt() > BEGIN) then

      call Flow % Grad_Variable(uu)
      call Flow % Grad_Variable(vv)
      call Flow % Grad_Variable(ww)
      call Flow % Grad_Variable(uv)
      call Flow % Grad_Variable(uw)
      call Flow % Grad_Variable(vw)
      call Flow % Grad_Variable(eps)
      if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
        call Flow % Grad_Variable(f22)
      end if

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            call Grid % Face_Normal(s, nx, ny, nz)

            phi_n = uu % x(c1) * nx + uu % y(c1) * ny + uu % z(c1) * nz
            uu % n(c2) = uu % n(c2) - bulk_vel * phi_n * dt

            phi_n = vv % x(c1) * nx + vv % y(c1) * ny + vv % z(c1) * nz
            vv % n(c2) = vv % n(c2) - bulk_vel * phi_n * dt

            phi_n = ww % x(c1) * nx + ww % y(c1) * ny + ww % z(c1) * nz
            ww % n(c2) = ww % n(c2) - bulk_vel * phi_n * dt

            phi_n = uv % x(c1) * nx + uv % y(c1) * ny + uv % z(c1) * nz
            uv % n(c2) = uv % n(c2) - bulk_vel * phi_n * dt

            phi_n = uw % x(c1) * nx + uw % y(c1) * ny + uw % z(c1) * nz
            uw % n(c2) = uw % n(c2) - bulk_vel * phi_n * dt

            phi_n = vw % x(c1) * nx + vw % y(c1) * ny + vw % z(c1) * nz
            vw % n(c2) = vw % n(c2) - bulk_vel * phi_n * dt

            phi_n = eps % x(c1) * nx + eps % y(c1) * ny + eps % z(c1) * nz
            eps % n(c2) = eps % n(c2) - bulk_vel * phi_n * dt

            if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
              phi_n = f22 % x(c1) * nx + f22 % y(c1) * ny + f22 % z(c1) * nz
              f22 % n(c2) = f22 % n(c2) - bulk_vel * phi_n * dt
            end if
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    else      ! Time % Curr_Dt() <= BEGIN

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            uu % n(c2) = uu % n(c1)
            vv % n(c2) = vv % n(c1)
            ww % n(c2) = ww % n(c1)
            uv % n(c2) = uv % n(c1)
            uw % n(c2) = uw % n(c1)
            vw % n(c2) = vw % n(c1)
            eps % n(c2) = eps % n(c1)
            if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
              f22 % n(c2) = f22 % n(c1)
            end if
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    end if    ! Time % Curr_Dt() > BEGIN

  end if

  !-------------!
  !             !
  !   Scalars   !
  !             !
  !-------------!

  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)

    if(Time % Curr_Dt() > BEGIN) then

      call Flow % Grad_Variable(phi)

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            call Grid % Face_Normal(s, nx, ny, nz)

            phi_n = phi % x(c1) * nx + phi % y(c1) * ny + phi % z(c1) * nz
            phi % n(c2) = phi % n(c2) - bulk_vel * phi_n * dt
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    else      ! Time % Curr_Dt() <= BEGIN

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            phi % n(c2) = phi % n(c1)
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    end if    ! Time % Curr_Dt() > BEGIN

  end do      ! sc

  !-------------------!
  !                   !
  !   Heat transfer   !
  !                   !
  !-------------------!

  if(Flow % heat_transfer) then

    if(Time % Curr_Dt() > BEGIN) then

      ! Temperature gradients might have been computed and
      ! stored already in t % x, t % y and t % z, check it
      call Flow % Grad_Variable(t)

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            call Grid % Face_Normal(s, nx, ny, nz)

            phi_n = t % x(c1) * nx + t % y(c1) * ny + t % z(c1) * nz
            t % n(c2) = t % n(c2) - bulk_vel * phi_n * dt
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    else      ! Time % Curr_Dt() <= BEGIN

      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. CONVECT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            t % n(c2) = t % n(c1)
          end do  ! face
        end if    ! boundary condition
      end do      ! region

    end if        ! Time % Curr_Dt() < BEGIN

  end if

  end subroutine
