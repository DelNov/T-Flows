!==============================================================================!
  subroutine Convective_Outflow(flow, turb, mult, dt)
!------------------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod
  use Turb_Mod
  use Grid_Mod
  use Control_Mod
  use Multiphase_Mod, only: Multiphase_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  real                          :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  type(Var_Type),  pointer :: kin, eps, zeta, f22, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Face_Type), pointer :: m_flux
  integer                  :: c1, c2, s
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  m_flux => flow % m_flux
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Field_Mod_Alias_Energy     (flow, t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_T2          (turb, t2)

  call Field_Mod_Calculate_Fluxes(flow, m_flux % n)

  call Field_Mod_Grad_Variable(flow, u)
  call Field_Mod_Grad_Variable(flow, v)
  call Field_Mod_Grad_Variable(flow, w)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! On the boundary perform the extrapolation
    if(c2  < 0) then
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
        u % n(c2) = u % n(c2)   &
                  - ( bulk % u * u % x(c1)         &
                    + bulk % v * u % y(c1)         &
                    + bulk % w * u % z(c1) ) * dt
        v % n(c2) = v % n(c2)  &
                  - ( bulk % u * v % x(c1)         &
                    + bulk % v * v % y(c1)         &
                    + bulk % w * v % z(c1) ) * dt
        w % n(c2) = w % n(c2)  &
                  - ( bulk % u * w % x(c1)         &
                    + bulk % v * w % y(c1)         &
                    + bulk % w * w % z(c1) ) * dt
      end if
    end if
  end do

  ! k-epsilon-zeta-f
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then

    call Field_Mod_Grad_Variable(flow, kin)
    call Field_Mod_Grad_Variable(flow, eps)
    call Field_Mod_Grad_Variable(flow, f22)
    call Field_Mod_Grad_Variable(flow, zeta)
    if(heat_transfer) then
      call Field_Mod_Grad_Variable(flow, t2)
    end if

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2  < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
          kin % n(c2) = kin % n(c2)   &
                      - ( bulk % u * kin % x(c1)         &
                        + bulk % v * kin % y(c1)         &
                        + bulk % w * kin % z(c1) ) * dt
          eps % n(c2) = eps % n(c2)   &
                      - ( bulk % u * eps % x(c1)         &
                        + bulk % v * eps % y(c1)         &
                        + bulk % w * eps % z(c1) ) * dt
          f22 % n(c2) = f22 % n(c2)   &
                      - ( bulk % u * f22 % x(c1)         &
                        + bulk % v * f22 % y(c1)         &
                        + bulk % w * f22 % z(c1) ) * dt
          zeta % n(c2) = zeta % n(c2)   &
                      - ( bulk % u * zeta % x(c1)         &
                        + bulk % v * zeta % y(c1)         &
                        + bulk % w * zeta % z(c1) ) * dt
          if(heat_transfer) then
            t2 % n(c2) = t2 % n(c2)   &
                        - ( bulk % u * t2 % x(c1)         &
                          + bulk % v * t2 % y(c1)         &
                          + bulk % w * t2 % z(c1) ) * dt
          end if
        end if
      end if
    end do

  end if

  if(heat_transfer) then

    ! Temperature gradients might have been computed and
    ! stored already in t % x, t % y and t % z, check it
    call Field_Mod_Grad_Variable(flow, t)

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2  < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
          t % n(c2) = t % n(c2)   &
                    - ( bulk % u * t % x(c1)         &
                      + bulk % v * t % y(c1)         &
                      + bulk % w * t % z(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
