!==============================================================================!
  subroutine Update_Boundary_Values(Flow, turb, Vof, update)
!------------------------------------------------------------------------------!
!   Update variables on the boundaries (boundary cells) where needed.          !
!   It does not update volume fluxes at boundaries - keep it like that.        !
!   Fluxes inside are calculated inside Compute_Pressure, and fluxes at        !
!   boundaries are updated in Balance_Volume                                   !
!                                                                              !
!   Ideally, this function should be called before each calculation of         !
!   gradients, but it is hardly so.  Maybe one could even think of calling     !
!   it from calculation of variables and updating those which match name       !
!   of the variable.                                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Field_Mod,      only: Field_Type
  use Turb_Mod
  use Grid_Mod
  use Control_Mod
  use Vof_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  character(*)             :: update
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, phi, fun
  type(Var_Type),  pointer :: kin, eps, zeta, f22, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  integer                  :: c1, c2, s, sc
  real                     :: kin_vis, u_tau
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  vis  => turb % vis
  fun  => Vof % fun
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_T2          (turb, t2)

  call To_Upper_Case(update)

  if(update .ne. 'MOMENTUM'   .and.  &
     update .ne. 'TURBULENCE' .and.  &
     update .ne. 'MULTIPHASE' .and.  &
     update .ne. 'ENERGY'     .and.  &
     update .ne. 'SCALARS'    .and.  &
     update .ne. 'ALL') then
    if(this_proc < 2) then
      print *, '# Invalid parameter in call to Update_Boundary_Values'
    end if
    call Comm_Mod_End
    stop
  end if

  !--------------!
  !              !
  !   Momentum   !
  !              !
  !--------------!
  if( (update .eq. 'MOMENTUM' .or. update .eq. 'ALL') ) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then

        ! Extrapolate velocities on the outflow boundary
        ! SYMMETRY is intentionally not treated here because I wanted to
        ! be sure that is handled only via graPHI and NewUVW functions)
        if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.     &
            Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.    &
            Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY ) then
          u % n(c2) = u % n(c1)
          v % n(c2) = v % n(c1)
          w % n(c2) = w % n(c1)
        end if
      end if ! c2 < 0
    end do

  end if  ! update momentum

  !----------------!
  !                !
  !   Multiphase   !
  !                !
  !----------------!
  if( (update .eq. 'MULTIPHASE' .or. update .eq. 'ALL') .and.  &
      Vof % model .eq. VOLUME_OF_FLUID ) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then

        if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.    &
            Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.    &
            Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.    &
            Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY ) then
          fun % n(c2) = fun % n(c1)
        end if
      end if ! c2 < 0
    end do

  end if  ! update multiphase

  !----------------!
  !                !
  !   Turbulence   !
  !                !
  !----------------!
  if(update .eq. 'TURBULENCE' .or. update .eq. 'ALL') then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then


        ! Spalart Allmaras
        if(turb % model .eq. SPALART_ALLMARAS .or.  &
           turb % model .eq. DES_SPALART) then
          if ( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.  &
               Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.  &
               Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.  &
               Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
            vis % n(c2) = vis % n(c1)
          end if
        end if

        ! Reynolds stress models
        if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
           turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
            uu  % n(c2) = 0.0
            vv  % n(c2) = 0.0
            ww  % n(c2) = 0.0
            uv  % n(c2) = 0.0
            uw  % n(c2) = 0.0
            vw  % n(c2) = 0.0
            kin % n(c2) = 0.0
            kin_vis = Flow % viscosity(c1) / Flow % density(c1)
            u_tau = kin_vis * sqrt(u % n(c1)**2 + v % n(c1)**2 + w % n(c1)**2)  &
                  / grid % wall_dist(c1)
            turb % y_plus(c1) = Y_Plus_Low_Re(turb, u_tau,           &
                                              grid % wall_dist(c1),  &
                                              kin_vis)
            if(turb % model .eq. RSM_MANCEAU_HANJALIC) f22 % n(c2) = 0.0
          end if
        end if

        ! k-epsilon-zeta-f
        if(turb % model .eq. K_EPS_ZETA_F .or.  &
           turb % model .eq. HYBRID_LES_RANS) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.   &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.   &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.   &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
            kin  % n(c2) = kin  % n(c1)
            eps  % n(c2) = eps  % n(c1)
            zeta % n(c2) = zeta % n(c1)
            f22  % n(c2) = f22  % n(c1)
            if(Flow % heat_transfer) then
              t2  % n(c2) = t2  % n(c1)
            end if
          end if

        end if

        ! k-epsilon
        if(turb % model .eq. K_EPS) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW  .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
            kin % n(c2) = kin % n(c1)
            eps % n(c2) = eps % n(c1)
            if(Flow % heat_transfer) then
              t2  % n(c2) = t2  % n(c1)
            end if 
          end if
        end if

        if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
           turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
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
            if(turb % model .eq. RSM_MANCEAU_HANJALIC)  &
              f22 % n(c2) = f22 % n(c1)
          end if
        end if
      end if ! c2 < 0
    end do

  end if  ! update turbulence

  !-------------------!
  !                   !
  !   Heat transfer   !
  !                   !
  !-------------------!
  if( (update .eq. 'ENERGY' .or. update .eq. 'ALL') .and.  &
      Flow % heat_transfer ) then

    ! Initialize variables for global heat transfer
    Flow % heat        = 0.0     ! [W]
    Flow % heat_flux   = 0.0     ! [W/m^2]
    Flow % heated_area = 0.0     ! [m^2]

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if (c2 < 0) then

        ! Wall temperature or heat fluxes for k-eps-zeta-f
        ! and high-re k-eps models. 
        if(turb % model .eq. K_EPS_ZETA_F    .or.  &
           turb % model .eq. HYBRID_LES_RANS .or.  &
           turb % model .eq. LES_DYNAMIC     .or.  &
           turb % model .eq. K_EPS) then
          if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * grid % wall_dist(c1)  &
                      / (turb % con_w(c1) + TINY)
          else if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * turb % con_w(c1)  &
                      / grid % wall_dist(c1)
          end if

        ! Wall temperature or heat fluxes for other trubulence models
        else
          if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * grid % wall_dist(c1)  &
                      / Flow % conductivity(c1)
          else if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * Flow % conductivity(c1)  &
                      / grid % wall_dist(c1)
          end if

        end if

        ! Integrate heat and heated area
        Flow % heat = Flow % heat + t % q(c2) * grid % s(s)
        if(abs(t % q(c2)) > TINY) then
          Flow % heated_area = Flow % heated_area + grid % s(s)
        end if

        if( Var_Mod_Bnd_Cond_Type(t,c2) .eq. OUTFLOW .or.     &
            Var_Mod_Bnd_Cond_Type(t,c2) .eq. PRESSURE .or.    &
            Var_Mod_Bnd_Cond_Type(t,c2) .eq. SYMMETRY ) then
          t % n(c2) = t % n(c1)
        end if

      end if ! c2 < 0
    end do ! s = 1, grid % n_faces

    !-----------------------------------------------!
    !   Integrate (summ) heated area, and heat up   !
    !-----------------------------------------------!
    call Comm_Mod_Global_Sum_Real(Flow % heat)
    call Comm_Mod_Global_Sum_Real(Flow % heated_area)
    Flow % heat_flux = Flow % heat / max(Flow % heated_area, TINY)

  end if  ! update energy and heat transfer

  !-------------!
  !             !
  !   Scalars   !
  !             !
  !-------------!
  if( (update .eq. 'SCALARS' .or. update .eq. 'ALL') ) then

    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)

      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        ! On the boundary perform the extrapolation
        if (c2 < 0) then

          ! Wall temperature or heat fluxes for k-eps-zeta-f
          ! and high-re k-eps models. 
          if(turb % model .eq. K_EPS_ZETA_F    .or.  &
             turb % model .eq. HYBRID_LES_RANS .or.  &
             turb % model .eq. K_EPS) then

            if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
              phi % n(c2) = phi % n(c1) + phi % q(c2) * grid % wall_dist(c1)  &
                        / (turb % diff_w(c1) + TINY)
            else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL) then
              phi % q(c2) = ( phi % n(c2) - phi % n(c1) ) * turb % diff_w(c1)  &
                        / grid % wall_dist(c1)
            end if

          ! Scalar boundary for other trubulence models
          else
            if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
              phi % n(c2) = phi % n(c1) + phi % q(c2) * grid % wall_dist(c1)  &
                          / Flow % diffusivity
            else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL) then
              phi % q(c2) = (phi % n(c2) - phi % n(c1)) * Flow % diffusivity &
                          / grid % wall_dist(c1)
            end if ! WALL or WALLFL
          end if ! Turb. models

          if( Var_Mod_Bnd_Cond_Type(phi,c2) .eq. OUTFLOW .or.     &
              Var_Mod_Bnd_Cond_Type(phi,c2) .eq. PRESSURE .or.    &
              Var_Mod_Bnd_Cond_Type(phi,c2) .eq. SYMMETRY ) then
            phi % n(c2) = phi % n(c1)
          end if

        end if ! c2 < 0
      end do ! s = 1, grid % n_faces
    end do ! sc = 1, Flow % n_scalars

  end if  ! update_scalars

  end subroutine
