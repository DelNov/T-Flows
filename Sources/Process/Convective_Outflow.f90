!==============================================================================!
  subroutine Convective_Outflow(flow, turb, Vof)
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
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t, phi
  type(Var_Type),  pointer :: kin, eps, zeta, f22, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Face_Type), pointer :: v_flux
  integer                  :: c1, c2, s, sc
  real                     :: dt
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  v_flux => flow % v_flux
  dt     =  flow % dt
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Field_Mod_Alias_Energy     (flow, t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_T2          (turb, t2)

  call Field_Mod_Calculate_Mass_Fluxes(flow, v_flux % n)

  !----------------------------------------------------------!
  !   Compute velocity gradients, taking jump into account   !
  !----------------------------------------------------------!

  !   ! If there is a jump in velocities, call specialized gradient calculation
  !   if(Vof % model .eq. VOLUME_OF_FLUID .and. flow % mass_transfer) then
  !     call Multiphase_Mod_Vof_Grad_Variable_With_Jump(Vof, flow % u)
  !     call Multiphase_Mod_Vof_Grad_Variable_With_Jump(Vof, flow % v)
  !     call Multiphase_Mod_Vof_Grad_Variable_With_Jump(Vof, flow % w)
  ! 
  !   ! No jumps, call usual routines
  !   else
  !     call Field_Mod_Grad_Variable(flow, flow % u)
  !     call Field_Mod_Grad_Variable(flow, flow % v)
  !     call Field_Mod_Grad_Variable(flow, flow % w)
  !   end if
  call Field_Mod_Grad_Variable(flow, flow % u)
  call Field_Mod_Grad_Variable(flow, flow % v)
  call Field_Mod_Grad_Variable(flow, flow % w)

  !------------------------!
  !                        !
  !   Momentum equations   !
  !                        !
  !------------------------!

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! On the boundary perform the extrapolation
    if(c2 < 0) then
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
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
  if(turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then

    call Field_Mod_Grad_Variable(flow, kin)
    call Field_Mod_Grad_Variable(flow, eps)
    if(flow % heat_transfer) then
      call Field_Mod_Grad_Variable(flow, t2)
    end if

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
          kin % n(c2) = kin % n(c2)                       &
                      - ( bulk % u * kin % x(c1)          &
                        + bulk % v * kin % y(c1)          &
                        + bulk % w * kin % z(c1) ) * dt
          eps % n(c2) = eps % n(c2)                       &
                      - ( bulk % u * eps % x(c1)          &
                        + bulk % v * eps % y(c1)          &
                        + bulk % w * eps % z(c1) ) * dt
          if(flow % heat_transfer) then
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

    call Field_Mod_Grad_Variable(flow, kin)
    call Field_Mod_Grad_Variable(flow, eps)
    call Field_Mod_Grad_Variable(flow, f22)
    call Field_Mod_Grad_Variable(flow, zeta)
    if(flow % heat_transfer) then
      call Field_Mod_Grad_Variable(flow, t2)
    end if

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
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
          if(flow % heat_transfer) then
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

    call Field_Mod_Grad_Variable(flow, uu)
    call Field_Mod_Grad_Variable(flow, vv)
    call Field_Mod_Grad_Variable(flow, ww)
    call Field_Mod_Grad_Variable(flow, uv)
    call Field_Mod_Grad_Variable(flow, uw)
    call Field_Mod_Grad_Variable(flow, vw)
    call Field_Mod_Grad_Variable(flow, eps)
    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Field_Mod_Grad_Variable(flow, f22)
    end if

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
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
  do sc = 1, flow % n_scalars
    phi => flow % scalar(sc)

    call Field_Mod_Grad_Variable(flow, phi)

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
          phi % n(c2) = phi % n(c2)                      &
                      - ( bulk % u * phi % x(c1)         &
                        + bulk % v * phi % y(c1)         &
                        + bulk % w * phi % z(c1) ) * dt
        end if
      end if
    end do
  end do

  if(flow % heat_transfer) then

    ! Temperature gradients might have been computed and
    ! stored already in t % x, t % y and t % z, check it
    call Field_Mod_Grad_Variable(flow, t)

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
          t % n(c2) = t % n(c2)                      &
                    - ( bulk % u * t % x(c1)         &
                      + bulk % v * t % y(c1)         &
                      + bulk % w * t % z(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
