!==============================================================================!
  subroutine Turb_Mod_Src_Rsm_Hanjalic_Jakirlic(turb, sol,  &
                                                name_phi, n_time_step)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for Hanjalic-Jakirlic model.                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: l_sc_x => r_cell_01,  &
                      l_sc_y => r_cell_02,  &
                      l_sc_z => r_cell_03,  &
                      kin_x  => r_cell_04,  &
                      kin_y  => r_cell_05,  &
                      kin_z  => r_cell_06,  &
                      kin_xx => r_cell_07,  &
                      kin_yy => r_cell_08,  &
                      kin_zz => r_cell_09,  &
                      ui_xx  => r_cell_10,  &
                      ui_yy  => r_cell_11,  &
                      ui_zz  => r_cell_12,  &
                      ui_xy  => r_cell_13,  &
                      ui_xz  => r_cell_14,  &
                      ui_yz  => r_cell_15,  &
                      kin_e  => r_cell_16,  &
                      kin_e_x=> r_cell_17,  &
                      kin_e_y=> r_cell_18,  &
                      kin_e_z=> r_cell_19,  &
                      diss1  => r_cell_20
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
  character(len=*)          :: name_phi
  integer                   :: n_time_step
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, zeta, f22
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, s, c1, c2, i, icont, nc, nb
  real                       :: mag
  real                       :: a11, a22, a33, a12, a13, a23
  real                       :: s11, s22, s33, s12, s13, s23
  real                       :: v11, v22, v33, v12, v13, v23
  real                       :: n1, n2, n3, aa2, aa3, aa, re_t
  real                       :: ff2, fd, ff1, cc, c1w, c2w, f_w, uu_nn
  real                       :: e11, e12, e13, e21, e22, e23, e31, e32, e33
  real                       :: eps11, eps12, eps13
  real                       :: eps21, eps22, eps23
  real                       :: eps31, eps32, eps33
  real                       :: f_eps, phi2_nn
  real                       :: fss,e2,e3,ee,cc1,cc2
  real                       :: uxx, uyy, uzz, uxy, uxz, uyz, uyx, uzx, uzy
  real                       :: r13, r23
  real                       :: a_lk_s_lk, a_mn_a_mn
  real                       :: var1w_11, var1w_22, var1w_33
  real                       :: var1w_12, var1w_13, var1w_23
  real                       :: var2w_11, var2w_22, var2w_33
  real                       :: var2w_12, var2w_13, var2w_23
  real                       :: var1_11,var1_22,var1_33,var1_12,var1_13,var1_23
  real                       :: var2_11,var2_22,var2_33,var2_12,var2_13,var2_23
  real                       :: p11, p22, p33, p12, p13, p23, eps_1, eps_2
  real                       :: ff5, tkolm, kin_vis
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   wall visc.    vis_wall [kg/(m*s)]  |                                       !
!   thermal cap.  capacity[m^2/(s^2*K)]| therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
! but dens > 1 mod. not applied here yet

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  nc   =  grid % n_cells
  nb   =  grid % n_bnd_cells
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Solver_Mod_Alias_System    (sol, a, b)

  diss1 = 0.0

  ee = 0.5
  aa = 0.5

  do c = 1, grid % n_cells
    kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
  end do

! !---------------------------------------------------!
! !   Below is one of versions of Hanjalic-Jakirlic   !
! !      model that required much more memory         !
! !---------------------------------------------------!
! if(name_phi == "23") then
!   call Field_Mod_Grad_Component(flow, uu % n, 1, var3x) ! duu/dx  
!   call Field_Mod_Grad_Component(flow, uu % n, 2, var3y) ! duu/dy  
!   call Field_Mod_Grad_Component(flow, uu % n, 3, var3z) ! duu/dz  
!
!   call Field_Mod_Grad_Component(flow, vv % n, 1, var4x) ! duw/dx  
!   call Field_Mod_Grad_Component(flow, vv % n, 2, var4y) ! duw/dy  
!   call Field_Mod_Grad_Component(flow, vv % n, 3, var4z) ! duw/dz  
!
!   call Field_Mod_Grad_Component(flow, ww % n, 1, var5x) ! duw/dx  
!   call Field_Mod_Grad_Component(flow, ww % n, 2, var5y) ! duw/dy  
!   call Field_Mod_Grad_Component(flow, ww % n, 3, var5z) ! duw/dz  
!
!   call Field_Mod_Grad_Component(flow, uv % n, 1, var6x) ! duv/dx  
!   call Field_Mod_Grad_Component(flow, uv % n, 2, var6y) ! duv/dy  
!   call Field_Mod_Grad_Component(flow, uv % n, 3, var6z) ! duv/dz  
!
!   call Field_Mod_Grad_Component(flow, uw % n, 1, kin_x) ! duw/dx  
!   call Field_Mod_Grad_Component(flow, uw % n, 2, kin_y) ! duw/dy  
!   call Field_Mod_Grad_Component(flow, uw % n, 3, kin_z) ! duw/dz  
!
!   call Field_Mod_Grad_Component(flow, vw % n, 1, var8x) ! duw/dx  
!   call Field_Mod_Grad_Component(flow, vw % n, 2, var8y) ! duw/dy  
!   call Field_Mod_Grad_Component(flow, vw % n, 3, var8z) ! duw/dz  
!
!   call Field_Mod_Grad_Component(flow, u % x, 1, var1x)  ! d2U/dxdx
!   call Field_Mod_Grad_Component(flow, u % y, 2, var1y)  ! d2U/dydy
!   call Field_Mod_Grad_Component(flow, u % z, 3, var1z)  ! d2U/dzdz
!   call Field_Mod_Grad_Component(flow, u % x, 2, var2x)  ! d2U/dxdy
!   call Field_Mod_Grad_Component(flow, u % x, 3, var2y)  ! d2U/dxdz
!   call Field_Mod_Grad_Component(flow, u % y, 3, var2z)  ! d2U/dydz
!
!   call Field_Mod_Grad_Component(flow, v % x, 1, var9x)  ! d2V/dxdx
!   call Field_Mod_Grad_Component(flow, v % y, 2, var9y)  ! d2V/dydy
!   call Field_Mod_Grad_Component(flow, v % z, 3, var9z)  ! d2V/dzdz
!   call Field_Mod_Grad_Component(flow, v % x, 2, var10x)  ! d2V/dxdy
!   call Field_Mod_Grad_Component(flow, v % x, 3, var10y)  ! d2V/dxdz
!   call Field_Mod_Grad_Component(flow, v % y, 3, var10z)  ! d2V/dydz
!
!   call Field_Mod_Grad_Component(flow, w % x, 1, var11x)  ! d2W/dxdx
!   call Field_Mod_Grad_Component(flow, w % y, 2, var11y)  ! d2W/dydy
!   call Field_Mod_Grad_Component(flow, w % z, 3, var11z)  ! d2W/dzdz
!   call Field_Mod_Grad_Component(flow, w % x, 2, var12x)  ! d2W/dxdy
!   call Field_Mod_Grad_Component(flow, w % x, 3, var12y)  ! d2W/dxdz
!   call Field_Mod_Grad_Component(flow, w % y, 3, var12z)  ! d2W/dydz
!
!   do c = 1, grid % n_cells
!     uxx = var1x(c)
!     uyy = var1y(c)
!     uzz = var1z(c)
!     uxy = var2x(c)
!     uxz = var2y(c)
!     uyz = var2z(c)
!     vxx = var9x(c)
!     vyy = var9y(c)
!     vzz = var9z(c)
!     vxy = var10x(c)
!     vxz = var10y(c)
!     Vyz = var10z(c)
!     wxx = var11x(c)
!     wyy = var11y(c)
!     wzz = var11z(c)
!     Wxy = var12x(c)
!     wxz = var12y(c)
!     wyz = var12z(c)
!     dudx= u % x(c) 
!     dudy= u % y(c) 
!     dudz= u % z(c) 
!     dvdx= v % x(c) 
!     dvdy= v % y(c) 
!     dvdz= v % z(c) 
!     dwdx= w % x(c) 
!     dwdy= w % y(c) 
!     dwdz= w % z(c) 
!     duu_dx = var3x(c)  
!     duu_dy = var3y(c)  
!     duu_dz = var3z(c)  
!     dvv_dx = var4x(c)  
!     dvv_dy = var4y(c)  
!     dvv_dz = var4z(c)  
!     dww_dx = var5x(c)  
!     dww_dy = var5y(c)  
!     dww_dz = var5z(c)  
!     duv_dx = var6x(c)  
!     duv_dy = var6y(c)  
!     duv_dz = var6z(c)  
!     duw_dx = kin_x(c)  
!     duw_dy = kin_y(c)  
!     duw_dz = kin_z(c)  
!     dvw_dx = var8x(c)  
!     dvw_dy = var8y(c)  
!     dvw_dz = var8z(c)  
!
!     diss1(c) = duu_dx*uxx + duv_dy*uyy + duw_dz*uzz   &
!              + uxy*(duv_dx + duu_dy)                  &
!              + uxz*(duw_dx + duu_dz)                  &
!              + uyz*(duw_dy + duv_dz)                  &
!              + duv_dx*vxx + dvv_dy*vyy + dvw_dz*vzz   &
!              + vxy*(dvv_dx + duv_dy)                  &
!              + vxz*(dvw_dx + duv_dz)                  &
!              + Vyz*(dvw_dy + dvv_dz)                  &
!              + duw_dx*wxx + dvw_dy*wyy + dww_dz*wzz   &
!              + Wxy*(dvw_dx + duw_dy)                  &
!              + wxz*(dww_dx + duw_dy)                  &
!              + wyz*(dww_dy + dvw_dz)                  &
!              + 0.32 * kin%n(c) / eps%n(c)  *  &
!                (uxx    * (duu_dx*dudx  + duv_dx*dudy  + duw_dx*dudz ) + & 
!                 uyy    * (duv_dy*dudx  + dvv_dy*dudy  + dvw_dy*dudz ) + & 
!                 uzz    * (duw_dz*dudx  + dvw_dz*dudy  + dww_dz*dudz ) + & 
!                 uxy    * (duu_dy*dudx  + duv_dy*dudy  + duw_dy*dudz   + &
!                           duv_dx*dudx  + dvv_dx*dudy  + dvw_dx*dudz ) + & 
!                 uxz    * (duu_dz*dudx  + duv_dz*dudy  + duw_dz*dudz   + &
!                           duw_dx*dudx  + dvw_dx*dudy  + dww_dx*dudz ) + &
!                 uyz    * (duv_dz*dudx  + dvv_dz*dudy  + dvw_dz*dudz   + &
!                           duw_dy*dudx  + dvw_dy*dudy  + dww_dy*dudz ) + &
!                 vxx    * (duu_dx*dvdx  + duv_dx*dvdy  + duw_dx*dvdz ) + & 
!                 vyy    * (duv_dy*dvdx  + dvv_dy*dvdy  + dvw_dy*dvdz ) + & 
!                 vzz    * (duw_dz*dvdx  + dvw_dz*dvdy  + dww_dz*dvdz ) + & 
!                 vxy    * (duu_dy*dvdx  + duv_dy*dvdy  + duw_dy*dvdz   + &
!                           duv_dx*dvdx  + dvv_dx*dvdy  + dvw_dx*dvdz ) + & 
!                 vxz    * (duu_dz*dvdx  + duv_dz*dvdy  + duw_dz*dvdz   + &
!                           duw_dx*dvdx  + dvw_dx*dvdy  + dww_dx*dvdz ) + &
!                 Vyz    * (duv_dz*dvdx  + dvv_dz*dvdy  + dvw_dz*dvdz   + &
!                           duw_dy*dvdx  + dvw_dy*dvdy  + dww_dy*dvdz ) + &
!                 wxx    * (duu_dx*dwdx  + duv_dx*dwdy  + duw_dx*dwdz ) + & 
!                 wyy    * (duv_dy*dwdx  + dvv_dy*dwdy  + dvw_dy*dwdz ) + & 
!                 wzz    * (duw_dz*dwdx  + dvw_dz*dwdy  + dww_dz*dwdz ) + & 
!                 Wxy    * (duu_dy*dwdx  + duv_dy*dwdy  + duw_dy*dwdz   + &
!                           duv_dx*dwdx  + dvv_dx*dwdy  + dvw_dx*dwdz ) + & 
!                 wxz    * (duu_dz*dwdx  + duv_dz*dwdy  + duw_dz*dwdz   + &
!                           duw_dx*dwdx  + dvw_dx*dwdy  + dww_dx*dwdz ) + &
!                 wyz    * (duv_dz*dwdx  + dvv_dz*dwdy  + dvw_dz*dwdz   + &
!                           duw_dy*dwdx  + dvw_dy*dwdy  + dww_dy*dwdz ))  
!     diss1(c) =  -2.0 * kin_vis * diss1(c)
!   end do
! end if

  if(name_phi == 'EPS') then
    do i=1,3
      if(i == 1) then
        call Field_Mod_Grad_Component(flow, u % x, 1, ui_xx(-nb:nc))  ! d2u/dxdx
        call Field_Mod_Grad_Component(flow, u % x, 2, ui_xy(-nb:nc))  ! d2u/dxdy
        call Field_Mod_Grad_Component(flow, u % x, 3, ui_xz(-nb:nc))  ! d2u/dxdz
        call Field_Mod_Grad_Component(flow, u % y, 2, ui_yy(-nb:nc))  ! d2u/dydy
        call Field_Mod_Grad_Component(flow, u % y, 3, ui_yz(-nb:nc))  ! d2u/dydz
        call Field_Mod_Grad_Component(flow, u % z, 3, ui_zz(-nb:nc))  ! d2u/dzdz
      end if
      if(i == 2) then
        call Field_Mod_Grad_Component(flow, v % x, 1, ui_xx(-nb:nc))  ! d2v/dxdx
        call Field_Mod_Grad_Component(flow, v % x, 2, ui_xy(-nb:nc))  ! d2v/dxdy
        call Field_Mod_Grad_Component(flow, v % x, 3, ui_xz(-nb:nc))  ! d2v/dxdz
        call Field_Mod_Grad_Component(flow, v % y, 2, ui_yy(-nb:nc))  ! d2v/dydy
        call Field_Mod_Grad_Component(flow, v % y, 3, ui_yz(-nb:nc))  ! d2v/dydz
        call Field_Mod_Grad_Component(flow, v % z, 3, ui_zz(-nb:nc))  ! d2v/dzdz
      end if
      if(i == 3) then
        call Field_Mod_Grad_Component(flow, w % x, 1, ui_xx(-nb:nc))  ! d2w/dxdx
        call Field_Mod_Grad_Component(flow, w % x, 2, ui_xy(-nb:nc))  ! d2w/dxdy
        call Field_Mod_Grad_Component(flow, w % x, 3, ui_xz(-nb:nc))  ! d2w/dxdz
        call Field_Mod_Grad_Component(flow, w % y, 2, ui_yy(-nb:nc))  ! d2w/dydy
        call Field_Mod_Grad_Component(flow, w % y, 3, ui_yz(-nb:nc))  ! d2w/dydz
        call Field_Mod_Grad_Component(flow, w % z, 3, ui_zz(-nb:nc))  ! d2w/dzdz
      end if

      do c = 1, grid % n_cells
        kin_vis = flow % viscosity(c) / flow % density(c)
        if(i == 1) then
          uxx = ui_xx(c)
          uxy = ui_xy(c)
          uyx = uxy
          uxz = ui_xz(c)
          uzx = uxz
          uyy = ui_yy(c)
          uyz = ui_yz(c)
          uzy = uyz
          uzz = ui_zz(c)
          diss1(c) =                                                           &
                  2.0 * 0.25 * kin_vis * kin % n(c) / max(eps % n(c), TINY)    &
                  * (  uu % n(c)*(uxx*uxx + uxy*uxy + uxz*uxz)                 &
                     + uv % n(c)*(uxx*uyx + uxy*uyy + uxz*uyz)                 &
                     + uw % n(c)*(uxx*uzx + uxy*uzy + uxz*uzz)                 &
                     + uv % n(c)*(uyx*uxx + uyy*uxy + uyz*uxz)                 &
                     + vv % n(c)*(uyx*uyx + uyy*uyy + uyz*uyz)                 &
                     + vw % n(c)*(uyx*uzx + uyy*uzy + uyz*uzz)                 &
                     + uw % n(c)*(uzx*uxx + uzy*uxy + uzz*uxz)                 &
                     + vw % n(c)*(uzx*uyx + uzy*uyy + uzz*uyz)                 &
                     + ww % n(c)*(uzx*uzx + uzy*uzy + uzz*uzz))
        end if
        if(i == 2) then
          uxx = ui_xx(c)
          uxy = ui_xy(c)
          uyx = uxy
          uxz = ui_xz(c)
          uzx = uxz
          uyy = ui_yy(c)
          uyz = ui_yz(c)
          uzy = uyz
          uzz = ui_zz(c)
          diss1(c) = diss1(c)                                 &
               +  2.0*0.25*kin_vis * kin % n(c) / eps % n(c)  &
               * (  uu % n(c)*(uxx*uxx + uxy*uxy + uxz*uxz)   &
                  + uv % n(c)*(uxx*uyx + uxy*uyy + uxz*uyz)   &
                  + uw % n(c)*(uxx*uzx + uxy*uzy + uxz*uzz)   &
                  + uv % n(c)*(uyx*uxx + uyy*uxy + uyz*uxz)   &
                  + vv % n(c)*(uyx*uyx + uyy*uyy + uyz*uyz)   &
                  + vw % n(c)*(uyx*uzx + uyy*uzy + uyz*uzz)   &
                  + uw % n(c)*(uzx*uxx + uzy*uxy + uzz*uxz)   &
                  + vw % n(c)*(uzx*uyx + uzy*uyy + uzz*uyz)   &
                  + ww % n(c)*(uzx*uzx + uzy*uzy + uzz*uzz))
        end if
        if(i == 3) then
          uxx = ui_xx(c)
          uxy = ui_xy(c)
          uyx = uxy
          uxz = ui_xz(c)
          uzx = uxz
          uyy = ui_yy(c)
          uyz = ui_yz(c)
          uzy = uyz
          uzz = ui_zz(c)
          diss1(c) = diss1(c)                                 &
               +  2.0*0.25*kin_vis * kin % n(c) / eps % n(c)  &
               * (  uu % n(c)*(uxx*uxx + uxy*uxy + uxz*uxz)   &
                  + uv % n(c)*(uxx*uyx + uxy*uyy + uxz*uyz)   &
                  + uw % n(c)*(uxx*uzx + uxy*uzy + uxz*uzz)   &
                  + uv % n(c)*(uyx*uxx + uyy*uxy + uyz*uxz)   &
                  + vv % n(c)*(uyx*uyx + uyy*uyy + uyz*uyz)   &
                  + vw % n(c)*(uyx*uzx + uyy*uzy + uyz*uzz)   &
                  + uw % n(c)*(uzx*uxx + uzy*uxy + uzz*uxz)   &
                  + vw % n(c)*(uzx*uyx + uzy*uyy + uzz*uyz)   &
                  + ww % n(c)*(uzx*uzx + uzy*uzy + uzz*uzz))
        end if
      end do
    end do  ! i
  end if    ! end if EPS == yes

  call Field_Mod_Grad_Component(flow, turb % l_scale, 1, l_sc_x(-nb:nc))
  call Field_Mod_Grad_Component(flow, turb % l_scale, 2, l_sc_y(-nb:nc))
  call Field_Mod_Grad_Component(flow, turb % l_scale, 3, l_sc_z(-nb:nc))

  r13 = ONE_THIRD
  r23 = TWO_THIRDS
  do  c = 1, grid % n_cells
    kin_vis = flow % viscosity(c) / flow % density(c)
    turb % p_kin(c) = max(                                                     &
          - (  uu % n(c)*u % x(c) + uv % n(c)*u % y(c) + uw % n(c)*u % z(c)    &
             + uv % n(c)*v % x(c) + vv % n(c)*v % y(c) + vw % n(c)*v % z(c)    &
             + uw % n(c)*w % x(c) + vw % n(c)*w % y(c) + ww % n(c)*w % z(c)),  &
               1.0e-10)

    mag = max(0.0, sqrt(l_sc_x(c)**2 + l_sc_y(c)**2 + l_sc_z(c)**2), TINY)

    n1 = l_sc_x(c) / mag
    n2 = l_sc_y(c) / mag
    n3 = l_sc_z(c) / mag

    a11 = uu % n(c) / kin % n(c) - r23
    a22 = vv % n(c) / kin % n(c) - r23
    a33 = ww % n(c) / kin % n(c) - r23
    a12 = uv % n(c) / kin % n(c)
    a13 = uw % n(c) / kin % n(c)
    a23 = vw % n(c) / kin % n(c)

    S11 = u % x(c)
    S22 = v % y(c)
    S33 = w % z(c)
    s12 = 0.5*(u % y(c)+v % x(c))
    s13 = 0.5*(u % z(c)+w % x(c))
    s23 = 0.5*(v % z(c)+w % y(c))

    V11 = 0.0
    V22 = 0.0
    V33 = 0.0
    V12 = 0.5*(u % y(c)-v % x(c)) - omega_z
    V13 = 0.5*(u % z(c)-w % x(c)) + omega_y
    V23 = 0.5*(v % z(c)-w % y(c)) - omega_x

    aa2 = (a11**2)+(a22**2)+(a33**2)+2*((a12**2)+(a13**2)+(a23**2))

    aa3 = a11**3 + a22**3 + a33**3 +                 &
          3*a12**2*(a11+a22) + 3*a13**2*(a11+a33) +  &
          3*a23**2*(a22+a33) + 6*a12*a13*a23

    aa=1.0 - (9.0/8.0)*(aa2-aa3)
    aa=max(aa,0.0)
    aa=min(aa,1.0)
 
    uu_nn = (uu % n(c)*n1*n1+uv % n(c)*n1*n2+uw % n(c)*n1*n3   &
           + uv % n(c)*n2*n1+vv % n(c)*n2*n2+vw % n(c)*n2*n3   &
           + uw % n(c)*n3*n1+vw % n(c)*n3*n2+ww % n(c)*n3*n3)

    a_mn_a_mn = a11*a11 + a22*a22 + a33*a33 + 2.0*(a12*a12+a13*a13+a23*a23)
    a_lk_s_lk = a11*S11 + a22*S22 + a33*S33 + 2.0*(a12*s12+a13*s13+a23*s23)

    re_t= (kin % n(c)**2)/(kin_vis*eps % n(c)+TINY)
    ff5 = min(aa2, (1.0-exp(-re_t/150))**3)
    tkolm = sqrt( kin_vis / max(eps % n(c), TINY) )

    ee=aa
    fss=1.0-(sqrt(aa) * ee**2)

    do icont=1,5
      eps11 = ((1.0 - fss)*r23 + fss*uu % n(c) / kin % n(c)  &
            + 2.0*ff5*S11*tkolm)*eps%n(c)
      eps22 = ((1.0 - fss)*r23 + fss*vv % n(c) / kin % n(c)  &
            + 2.0*ff5*S22*tkolm)*eps%n(c)
      eps33 = ((1.0 - fss)*r23 + fss*ww % n(c) / kin % n(c)  &
            + 2.0*ff5*S33*tkolm)*eps%n(c)
      eps12 = (fss*uv % n(c) / kin % n(c) + 2.0*ff5*s12*tkolm)*eps % n(c)
      eps13 = (fss*uw % n(c) / kin % n(c) + 2.0*ff5*s13*tkolm)*eps % n(c)
      eps23 = (fss*vw % n(c) / kin % n(c) + 2.0*ff5*s23*tkolm)*eps % n(c)
      eps21 = eps12
      eps31 = eps13
      eps32 = eps23

      e11 = eps11 / max(eps % n(c), TINY) - r23
      e22 = eps22 / max(eps % n(c), TINY) - r23
      e33 = eps33 / max(eps % n(c), TINY) - r23
      e12 = eps12 / max(eps % n(c), TINY)
      e13 = eps13 / max(eps % n(c), TINY)
      e23 = eps23 / max(eps % n(c), TINY)
      e21 = e12
      e31 = e13
      e32 = e23
      e2 = (e11**2)+(e22**2)+(e33**2)+2*((e12**2)+(e13**2)+(e23**2))

      e3 =   e11**3 + e22**3 + e33**3  &
         + 3*e12**2*(e11+e22)          &
         + 3*e13**2*(e11+e33)          &
         + 3*e23**2*(e22+e33)          &
         + 6*e12*e13*e23

      ee = 1.0 - (9.0/8.0)*(e2-e3)

      ee = max(ee, 0.0)
      ee = min(ee, 1.0)
      fss=1.0-(aa**0.5*ee**2)
    end do

    re_t  = (kin % n(c)*kin % n(c)) / (kin_vis*eps % n(c) + TINY)
    f_eps = 1.0 - ((c_2e-1.4)/c_2e)*exp(-(re_t/6.0)**2)
    ff2   = (re_t/150.0)**1.5
    ff2   = min(ff2,1.0)
    fd    = 1.0/(1.0+0.1*re_t)
    FF1   = min(0.6, aa2)
    cc    = 2.5*aa*FF1**0.25*ff2
    cc1   = cc+SQRT(aa)*(ee**2)
    cc2   = 0.8*SQRT(aa)
    c1w   = max((1.0 - 0.7*cc), 0.3)
    c2w   = min(aa,0.3)
    f_w   = min( kin % n(c)**1.5                                         &
                 / (2.5 * max(eps % n(c), TINY) * grid % wall_dist(c)),  &
                 1.4)

    p11 = - 2.0*(  uu % n(c) * u % x(c)      &
                 + uv % n(c) * u % y(c)      &
                 + uw % n(c) * u % z(c))     &
          - 2.0 * omega_y * 2.0 * uw % n(c)  &
          + 2.0 * omega_z * 2.0 * uv % n(c) 

    p22 = - 2.0*(  uv % n(c) * v % x(c)      &
                 + vv % n(c) * v % y(c)      &
                 + vw % n(c) * v % z(c))     &
          + 2.0 * omega_x * 2.0 * vw % n(c)  &
          - 2.0 * omega_z * 2.0 * uw % n(c) 

    p33 = - 2.0*(  uw % n(c) * w % x(c)      &
                 + vw % n(c) * w % y(c)      &
                 + ww % n(c) * w % z(c))     &
          - 2.0 * omega_x * 2.0 * vw % n(c)  &
          + 2.0 * omega_y * 2.0 * uw % n(c) 

    p12 = -(  uu % n(c) * v % x(c)       &
            + uv % n(c) * v % y(c)       &
            + uw % n(c) * v % z(c)       &
            + uv % n(c) * u % x(c)       &
            + vv % n(c) * u % y(c)       &
            + vw % n(c) * u % z(c))      &
            + 2.0 * omega_x * uw % n(c)  &
            - 2.0 * omega_y * vw % n(c)  &
            + 2.0 * omega_z * (vv % n(c) - uu % n(c)) 

    p13 = -(  uw % n(c)*u % x(c)                       &
            + vw % n(c)*u % y(c)                       &
            + ww % n(c)*u % z(c)                       &
            + uu % n(c)*w % x(c)                       &
            + uv % n(c)*w % y(c)                       &
            + uw % n(c)*w % z(c))                      &
            - 2.0 * omega_x * uv % n(c)                &
            - 2.0 * omega_y * (ww % n(c) - uu % n(c))  &
            + 2.0 * omega_z * vw % n(c) 

    p23 = -(  uv % n(c) * w % x(c)                     &
            + vv % n(c) * w % y(c)                     &
            + vw % n(c) * w % z(c)                     &
            + uw % n(c) * v % x(c)                     &
            + vw % n(c) * v % y(c)                     &
            + ww % n(c) * v % z(c))                    &
            - 2.0 * omega_x * (vw % n(c) - ww % n(c))  &
            + 2.0 * omega_y * uv % n(c)                &
            - 2.0 * omega_z * uw % n(c) 

    var1_11 = -cc1*eps%n(c)*a11 
    var1_22 = -cc1*eps%n(c)*a22 
    var1_33 = -cc1*eps%n(c)*a33 
    var1_12 = -cc1*eps%n(c)*a12 
    var1_13 = -cc1*eps%n(c)*a13 
    var1_23 = -cc1*eps%n(c)*a23 

    var2_11 = -cc2*(p11 - r23 * turb % p_kin(c))
    var2_22 = -cc2*(p22 - r23 * turb % p_kin(c))
    var2_33 = -cc2*(p33 - r23 * turb % p_kin(c))
    var2_12 = -cc2*p12
    var2_13 = -cc2*p13
    var2_23 = -cc2*p23

    phi2_nn =   var2_11 * n1 * n1  &
            + 2*var2_12 * n1 * n2  &
            + 2*var2_13 * n1 * n3  &
            +   var2_22 * n2 * n2  &
            + 2*var2_23 * n2 * n3  &
            +   var2_33 * n3 * n3

    var1w_11 = c1w * f_w * eps % n(c) / kin % n(c)     &
             * (uu_nn-1.5*2.0*(  uu % n(c)*n1*n1*0.0   &
                               + uv % n(c)*n1*n2       &
                               + uw % n(c)*n1*n3))
    var1w_22 = c1w * f_w * eps % n(c) / kin % n(c)     &
             * (uu_nn-1.5*2.0*(  uv % n(c)*n2*n1       &
                               + vv % n(c)*n2*n2*0.0   &
                               + vw % n(c)*n2*n3))
    var1w_33 = c1w * f_w * eps % n(c) / kin % n(c)     &
             * (uu_nn-1.5*2.0*(  uw % n(c)*n3*n1       &
                               + vw % n(c)*n3*n2       &
                               + ww % n(c)*n3*n3*0.0))
    var1w_12 = c1w * f_w * eps % n(c) / kin % n(c)  &
             * (-1.5*(  uu % n(c)*n2*n1             &
                      + uv % n(c)*n2*n2*0.0         &
                      + uw % n(c)*n2*n3             &
                      + uv % n(c)*n1*n1*0.0         &
                      + vv % n(c)*n1*n2             &
                      + vw % n(c)*n1*n3))
    var1w_13 = c1w * f_w * eps % n(c) / kin % n(c)  &
             *(-1.5*(  uu % n(c)*n3*n1              &
                     + uv % n(c)*n3*n2              &
                     + uw % n(c)*n3*n3*0.0          &
                     + uw % n(c)*n1*n1*0.0          &
                     + vw % n(c)*n1*n2              &
                     + ww % n(c)*n1*n3))
    var1w_23 = c1w * f_w * eps % n(c) / kin % n(c)  &
             *(-1.5*(  uw % n(c)*n2*n1              &
                     + vw % n(c)*n2*n2*0.0          &
                     + ww % n(c)*n2*n3              &
                     + uv % n(c)*n3*n1              &
                     + vv % n(c)*n3*n2              &
                     + vw % n(c)*n3*n3)*0.0)

    var2w_11 = c2w * f_w * (phi2_nn - 1.5*2.0*(  var2_11*n1*n1    &
                                               + var2_12*n1*n2    &
                                               + var2_13*n1*n3))
    var2w_22 = c2w * f_w * (phi2_nn - 1.5*2.0*(  var2_12*n1*n2    &
                                               + var2_22*n2*n2    &
                                               + var2_23*n3*n2))
    var2w_33 = c2w * f_w * (phi2_nn - 1.5*2.0*(  var2_13*n1*n3    &
                                               + var2_23*n2*n3    &
                                               + var2_33*n3*n3))
    var2w_12 = c2w * f_w * (-1.5*(  var2_11*n2*n1    &
                                  + var2_12*n2*n2    &
                                  + var2_13*n2*n3    &
                                  + var2_12*n1*n1    &
                                  + var2_22*n1*n2    &
                                  + var2_23*n1*n3))
    var2w_13 = c2w * f_w * (-1.5*(  var2_11*n3*n1    &
                                  + var2_12*n3*n2    &
                                  + var2_13*n3*n3    &
                                  + var2_13*n1*n1    &
                                  + var2_23*n1*n2    &
                                  + var2_33*n1*n3))
    var2w_23 = c2w * f_w * (-1.5*(  var2_13*n2*n1    &
                                  + var2_23*n2*n2    &
                                  + var2_33*n2*n3    &
                                  + var2_12*n3*n1    &
                                  + var2_22*n3*n2    &
                                  + var2_23*n3*n3))

    !---------------!
    !   UU stress   !
    !---------------!
    if(name_phi == 'UU') then

      b(c) = b(c) + flow % density(c) * (  max(p11,0.0)  &
                     + cc1 * eps % n(c) * r23            &
                     + max(var2_11, 0.0)                 &
                     + max(var1w_11,0.0)                 &
                     + max(var2w_11,0.0))*grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))                               &
                + flow % density(c) * (  cc1 * eps % n(c) / kin % n(c)        &
                             + c1w * f_w * eps % n(c) / kin % n(c)*3.0*n1*n1  &
                             + fss * eps % n(c) / kin % n(c))*grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))               &
                + flow % density(c) * (  max(-p11,     0.0)   &
                             + max(-var2_11, 0.0)             &
                             + max(-var1w_11,0.0)             &
                             + max(-var2w_11,0.0)             &
                             + (1.0-fss) * r23 * eps % n(c))  &
                            / max(uu%n(c), TINY) * grid % vol(c)

    !---------------!
    !   VV stress   !
    !---------------!
    else if(name_phi == 'VV') then

      b(c) = b(c) + flow % density(c) * (  max(p22,0.0)  &
                     + cc1 * eps % n(c) * r23            &
                     + max(var2_22, 0.0)                 &
                     + max(var1w_22,0.0)                 &
                     + max(var2w_22,0.0))*grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))                               &
                + flow % density(c) * (  cc1 * eps % n(c) / kin % n(c)        &
                             + c1w * f_w * eps % n(c) / kin % n(c)*3.0*n2*n2  &
                             + fss * eps % n(c) / kin % n(c))*grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))               &
                + flow % density(c) * (  max(-p22,     0.0)   &
                             + max(-var2_22, 0.0)             &
                             + max(-var1w_22,0.0)             &
                             + max(-var2w_22,0.0)             &
                             + (1.0-fss) * r23 * eps % n(c))  &
                            / max(vv%n(c), TINY) * grid % vol(c)

    !---------------!
    !   WW stress   !
    !---------------!
    else if(name_phi == 'WW') then

      b(c) = b(c) + flow % density(c) * (  max(p33,0.0)  &
                     + cc1 * eps % n(c) * r23            &
                     + max(var2_33, 0.0)                 &
                     + max(var1w_33,0.0)                 &
                     + max(var2w_33,0.0))*grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))                               &
                + flow % density(c) * (  cc1 * eps % n(c) / kin % n(c)        &
                             + c1w * f_w * eps % n(c) / kin % n(c)*3.0*n3*n3  &
                             + fss * eps % n(c) / kin % n(c))*grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))               &
                + flow % density(c) * (  max(-p33,     0.0)   &
                             + max(-var2_33, 0.0)             &
                             + max(-var1w_33,0.0)             &
                             + max(-var2w_33,0.0)             &
                             + (1.0-fss) * r23 * eps % n(c))  &
                            / max(ww % n(c), TINY) * grid % vol(c)

    !---------------!
    !   UV stress   !
    !---------------!
    else if(name_phi == 'UV') then
      b(c) = b(c) + flow % density(c) * (p12 + var2_12 + var1w_12 + var2w_12)  &
           * grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))                             &
            + flow % density(c) * (  cc1 * eps % n(c) / kin % n(c)          &
               + c1w * f_w * eps % n(c) / kin % n(c) * 1.5*(n1*n1 + n2*n2)  &
               + fss * eps % n(c) / kin % n(c) ) * grid % vol(c)

    !---------------!
    !   UW stress   !
    !---------------!
    else if(name_phi == 'UW') then
      b(c) = b(c) + flow % density(c) * (p13 + var2_13 + var1w_13 + var2w_13)  &
           * grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))                             &
            + flow % density(c) * (  cc1 * eps % n(c) / kin % n(c)          &
               + c1w * f_w * eps % n(c) / kin % n(c) * 1.5*(n1*n1 + n3*n3)  &
               + fss * eps % n(c) / kin % n(c) ) * grid % vol(c)

    !---------------!
    !   VW stress   !
    !---------------!
    else if(name_phi == 'VW') then
      b(c) = b(c) + flow % density(c) * (p23 + var2_23 + var1w_23 + var2w_23)  &
           * grid % vol(c)
      a % val(a % dia(c)) = a % val(a % dia(c))                             &
            + flow % density(c) * (  cc1 * eps % n(c) / kin % n(c)          &
               + c1w * f_w * eps % n(c) / kin % n(c) * 1.5*(n2*n2 + n3*n3)  &
               + fss * eps % n(c) / kin % n(c) ) * grid % vol(c)

    !----------------------!
    !   Epsilon equation   !
    !----------------------!
    else if(name_phi == 'EPS') then 
      f_eps = 1.0 - ((c_2e-1.4)/c_2e) * exp(-(re_t/6.0)**2)
      eps_1 = 1.44 * turb % p_kin(c) / turb % t_scale(c)
      eps_2 = c_2e * f_eps  / turb % t_scale(c)
      b(c) = b(c) + flow % density(c) * (eps_1 + diss1(c)) * grid % vol(c)

      a % val(a % dia(c)) = a % val(a % dia(c)) + flow % density(c) * eps_2  &
                          * grid % vol(c)
    end if
  end do

  do c = 1, grid % n_cells
    kin_e(c) = sqrt( 0.5 * (uu % n(c) + vv % n(c) + ww % n(c)) )
  end do

  if(name_phi == 'EPS') then
    call Field_Mod_Grad_Component(flow, kin_e(-nb:nc), 1, kin_e_x(-nb:nc))
    call Field_Mod_Grad_Component(flow, kin_e(-nb:nc), 2, kin_e_y(-nb:nc))
    call Field_Mod_Grad_Component(flow, kin_e(-nb:nc), 3, kin_e_z(-nb:nc))
    do c = 1, grid % n_cells
      kin_vis = flow % viscosity(c) / flow % density(c)
      re_t  = (kin % n(c)**2) / (kin_vis*eps % n(c) + TINY)
      f_eps = 1.0 - ((c_2e-1.4)/c_2e) * exp(-(re_t/6.0)**2)
      b(c) = b(c) + flow % density(c)                               &
                     * (c_2e * f_eps / turb % t_scale(c)            &
                     * (kin_vis *(  kin_e_x(c)**2                   &
                                  + kin_e_y(c)**2                   &
                                  + kin_e_z(c)**2)))                &
                     * grid % vol(c)
    end do
  end if

  if(name_phi == 'EPS') then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Calculate a values of dissipation  on wall
      if(c2 < 0 ) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

          kin_vis = flow % viscosity(c1) / flow % density(c1)
          eps % n(c2) = kin_vis * (  kin_e_x(c1)**2  &
                                   + kin_e_y(c1)**2  &
                                   + kin_e_z(c1)**2)

        end if   ! end if of BC=wall
      end if    ! end if of c2<0
    end do
  end if

  end subroutine
