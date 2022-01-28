!==============================================================================!
  subroutine Convective_Outflow(Flow, turb, Vof)
!------------------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod
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
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t, phi
  type(Var_Type),  pointer :: kin, eps, zeta, f22, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Face_Type), pointer :: v_flux
  integer                  :: c1, c2, s, sc
  real                     :: dt
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  v_flux => Flow % v_flux
  dt     =  Flow % dt
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_T2          (turb, t2)

  call Flow % Calculate_Fluxes(v_flux % n)

  !--------------------------------!
  !   Compute velocity gradients   !
  !--------------------------------!
  call Flow % Grad_Variable(Flow % u)
  call Flow % Grad_Variable(Flow % v)
  call Flow % Grad_Variable(Flow % w)

  !------------------------!
  !                        !
  !   Momentum equations   !
  !                        !
  !------------------------!

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! On the boundary perform the extrapolation
    if(c2 < 0) then
      if( (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) ) then
        u % n(c2) = u % n(c2)                      &
                  - ( bulk % u * u % x(c1)         &
                    + bulk % v * u % y(c1)         &
                    + bulk % w * u % z(c1) ) * dt
        v % n(c2) = v % n(c2)                      &
                  - ( bulk % u * v % x(c1)         &
                    + bulk % v * v % y(c1)         &
                    + bulk % w * v % z(c1) ) * dt
        w % n(c2) = w % n(c2)                      &
                  - ( bulk % u * w % x(c1)         &
                    + bulk % v * w % y(c1)         &
                    + bulk % w * w % z(c1) ) * dt
      end if
    end if
  end do

  !--------------------------!
  !                          !
  !   Turbulent quantities   !
  !                          !
  !--------------------------!

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(turb % model .eq. K_EPS) then

    call Flow % Grad_Variable(kin)
    call Flow % Grad_Variable(eps)
    if(Flow % heat_transfer) then
      call Flow % Grad_Variable(t2)
    end if

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) ) then
          kin % n(c2) = kin % n(c2)                       &
                      - ( bulk % u * kin % x(c1)          &
                        + bulk % v * kin % y(c1)          &
                        + bulk % w * kin % z(c1) ) * dt
          eps % n(c2) = eps % n(c2)                       &
                      - ( bulk % u * eps % x(c1)          &
                        + bulk % v * eps % y(c1)          &
                        + bulk % w * eps % z(c1) ) * dt
          if(Flow % heat_transfer) then
            t2 % n(c2) = t2 % n(c2)                       &
                        - ( bulk % u * t2 % x(c1)         &
                          + bulk % v * t2 % y(c1)         &
                          + bulk % w * t2 % z(c1) ) * dt
          end if
        end if
      end if
    end do

  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then

    call Flow % Grad_Variable(kin)
    call Flow % Grad_Variable(eps)
    call Flow % Grad_Variable(f22)
    call Flow % Grad_Variable(zeta)
    if(Flow % heat_transfer) then
      call Flow % Grad_Variable(t2)
    end if

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) ) then
          kin % n(c2) = kin % n(c2)                       &
                      - ( bulk % u * kin % x(c1)          &
                        + bulk % v * kin % y(c1)          &
                        + bulk % w * kin % z(c1) ) * dt
          eps % n(c2) = eps % n(c2)                       &
                      - ( bulk % u * eps % x(c1)          &
                        + bulk % v * eps % y(c1)          &
                        + bulk % w * eps % z(c1) ) * dt
          f22 % n(c2) = f22 % n(c2)                       &
                      - ( bulk % u * f22 % x(c1)          &
                        + bulk % v * f22 % y(c1)          &
                        + bulk % w * f22 % z(c1) ) * dt
          zeta % n(c2) = zeta % n(c2)                     &
                      - ( bulk % u * zeta % x(c1)         &
                        + bulk % v * zeta % y(c1)         &
                        + bulk % w * zeta % z(c1) ) * dt
          if(Flow % heat_transfer) then
            t2 % n(c2) = t2 % n(c2)                       &
                        - ( bulk % u * t2 % x(c1)         &
                          + bulk % v * t2 % y(c1)         &
                          + bulk % w * t2 % z(c1) ) * dt
          end if
        end if
      end if
    end do

  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!

  if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

    call Flow % Grad_Variable(uu)
    call Flow % Grad_Variable(vv)
    call Flow % Grad_Variable(ww)
    call Flow % Grad_Variable(uv)
    call Flow % Grad_Variable(uw)
    call Flow % Grad_Variable(vw)
    call Flow % Grad_Variable(eps)
    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Flow % Grad_Variable(f22)
    end if

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) ) then
          uu % n(c2) = uu % n(c2)                          &
                      - ( bulk % u * uu % x(c1)            &
                        + bulk % v * uu % y(c1)            &
                        + bulk % w * uu % z(c1) ) * dt
          vv % n(c2) = vv % n(c2)                          &
                      - ( bulk % u * vv % x(c1)            &
                        + bulk % v * vv % y(c1)            &
                        + bulk % w * vv % z(c1) ) * dt
          ww % n(c2) = ww % n(c2)                          &
                      - ( bulk % u * ww % x(c1)            &
                        + bulk % v * ww % y(c1)            &
                        + bulk % w * ww % z(c1) ) * dt
          uv % n(c2) = uv % n(c2)                          &
                      - ( bulk % u * uv % x(c1)            &
                        + bulk % v * uv % y(c1)            &
                        + bulk % w * uv % z(c1) ) * dt
          uw % n(c2) = uw % n(c2)                          &
                      - ( bulk % u * uw % x(c1)            &
                        + bulk % v * uw % y(c1)            &
                        + bulk % w * uw % z(c1) ) * dt
          vw % n(c2) = vw % n(c2)                          &
                      - ( bulk % u * vw % x(c1)            &
                        + bulk % v * vw % y(c1)            &
                        + bulk % w * vw % z(c1) ) * dt
          eps % n(c2) = eps % n(c2)                        &
                      - ( bulk % u * eps % x(c1)           &
                        + bulk % v * eps % y(c1)           &
                        + bulk % w * eps % z(c1) ) * dt
          if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
            f22 % n(c2) = f22 % n(c2)                      &
                        - ( bulk % u * f22 % x(c1)         &
                          + bulk % v * f22 % y(c1)         &
                          + bulk % w * f22 % z(c1) ) * dt
          end if
        end if
      end if
    end do

  end if

  !-------------!
  !   Scalars   !
  !-------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)

    call Flow % Grad_Variable(phi)

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) ) then
          phi % n(c2) = phi % n(c2)                      &
                      - ( bulk % u * phi % x(c1)         &
                        + bulk % v * phi % y(c1)         &
                        + bulk % w * phi % z(c1) ) * dt
        end if
      end if
    end do
  end do

  if(Flow % heat_transfer) then

    ! Temperature gradients might have been computed and
    ! stored already in t % x, t % y and t % z, check it
    call Flow % Grad_Variable(t)

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) ) then
          t % n(c2) = t % n(c2)                      &
                    - ( bulk % u * t % x(c1)         &
                      + bulk % v * t % y(c1)         &
                      + bulk % w * t % z(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
