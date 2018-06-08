!==============================================================================!
  subroutine Source_Hanjalic_Jakirlic(grid, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for Hanjalic-Jakirlic model.                                   !  
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
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
                      diss1  => r_cell_17    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2, i, icont
  real    :: mag
  real    :: a11, a22, a33, a12, a13, a23
  real    :: s11, s22, s33, s12, s13, s23
  real    :: v11, v22, v33, v12, v13, v23
  real    :: n1,n2,n3,aa2,aa3,aa,re_t,ff2,fd,ff1,CC,c1w,c2w,f_w,uu_nn
  real    :: e11,e12,e13,e21,e22,e23,e31,e32,e33
  real    :: eps11,eps12,eps13,eps21,eps22,eps23,eps31,eps32,eps33
  real    :: f_eps, phi2_nn, eps_over_kin
  real    :: fss,E2,E3,ee,cc1,cc2
  real    :: Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, Uyx, Uzx, Uzy
  real    :: r13, r23
  real    :: a_lk_s_lk, a_mn_a_mn
  real    :: var1w_11, var1w_22, var1w_33, var1w_12, var1w_13, var1w_23
  real    :: var2w_11, var2w_22, var2w_33, var2w_12, var2w_13, var2w_23
  real    :: var1_11, var1_22, var1_33, var1_12, var1_13, var1_23
  real    :: var2_11, var2_22, var2_33, var2_12, var2_13, var2_23
  real    :: p11, p22, p33, p12, p13, p23, eps_1, eps_2
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   Wall visc.    vis_wall [kg/(m*s)]  |                                       !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
! but dens > 1 mod. not applied here yet

  diss1 = 0.0

  ee = 0.5
  aa = 0.5

  do c = 1, grid % n_cells
    kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
    l_scale(c) = kin % n(c)**1.5/(eps % n(c)+TINY)
    t_scale(c) = kin % n(c)/(eps % n(c)+TINY)
  end do

  call Grad_Mod_For_Phi(grid, kin % n, 1, kin_x, .true.)  ! dK/dx
  call Grad_Mod_For_Phi(grid, kin % n, 2, kin_y, .true.)  ! dK/dy
  call Grad_Mod_For_Phi(grid, kin % n, 3, kin_z, .true.)  ! dK/dz

  call Grad_Mod_For_Phi(grid, kin_x, 1, kin_xx, .true.)  ! d^2 K / dx^2
  call Grad_Mod_For_Phi(grid, kin_y, 2, kin_yy, .true.)  ! d^2 K / dy^2
  call Grad_Mod_For_Phi(grid, kin_z, 3, kin_zz, .true.)  ! d^2 K / dz^2

  do c = 1, grid % n_cells
    eps_tot(c) = max(eps % n(c) + &
      0.5 * viscosity * (kin_xx(c) + kin_yy(c) + kin_zz(c)), TINY)
  end do

! !---------------------------------------------------!
! !   Below is one of versions of Hanjalic-Jakirlic   !
! !      model that required much more memory         !
! !---------------------------------------------------!
! if(name_phi .eq. "23") then
!   call Grad_Mod_For_Phi(grid, uu % n, 1, var3x, .true.) ! duu/dx  
!   call Grad_Mod_For_Phi(grid, uu % n, 2, var3y, .true.) ! duu/dy  
!   call Grad_Mod_For_Phi(grid, uu % n, 3, var3z, .true.) ! duu/dz  
!
!   call Grad_Mod_For_Phi(grid, vv % n, 1, var4x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, vv % n, 2, var4y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, vv % n, 3, var4z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, ww % n, 1, var5x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, ww % n, 2, var5y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, ww % n, 3, var5z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, uv % n, 1, var6x, .true.) ! duv/dx  
!   call Grad_Mod_For_Phi(grid, uv % n, 2, var6y, .true.) ! duv/dy  
!   call Grad_Mod_For_Phi(grid, uv % n, 3, var6z, .true.) ! duv/dz  
!
!   call Grad_Mod_For_Phi(grid, uw % n, 1, kin_x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, uw % n, 2, kin_y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, uw % n, 3, kin_z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, vw % n, 1, var8x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, vw % n, 2, var8y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, vw % n, 3, var8z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, u % x, 1, var1x, .true.)  ! d2U/dxdx
!   call Grad_Mod_For_Phi(grid, u % y, 2, var1y, .true.)  ! d2U/dydy
!   call Grad_Mod_For_Phi(grid, u % z, 3, var1z, .true.)  ! d2U/dzdz
!   call Grad_Mod_For_Phi(grid, u % x, 2, var2x, .true.)  ! d2U/dxdy
!   call Grad_Mod_For_Phi(grid, u % x, 3, var2y, .true.)  ! d2U/dxdz
!   call Grad_Mod_For_Phi(grid, u % y, 3, var2z, .true.)  ! d2U/dydz
!
!   call Grad_Mod_For_Phi(grid, v % x, 1, var9x, .true.)  ! d2V/dxdx
!   call Grad_Mod_For_Phi(grid, v % y, 2, var9y, .true.)  ! d2V/dydy
!   call Grad_Mod_For_Phi(grid, v % z, 3, var9z, .true.)  ! d2V/dzdz
!   call Grad_Mod_For_Phi(grid, v % x, 2, var10x, .true.)  ! d2V/dxdy
!   call Grad_Mod_For_Phi(grid, v % x, 3, var10y, .true.)  ! d2V/dxdz
!   call Grad_Mod_For_Phi(grid, v % y, 3, var10z, .true.)  ! d2V/dydz
!
!   call Grad_Mod_For_Phi(grid, w % x, 1, var11x, .true.)  ! d2W/dxdx
!   call Grad_Mod_For_Phi(grid, w % y, 2, var11y, .true.)  ! d2W/dydy
!   call Grad_Mod_For_Phi(grid, w % z, 3, var11z, .true.)  ! d2W/dzdz
!   call Grad_Mod_For_Phi(grid, w % x, 2, var12x, .true.)  ! d2W/dxdy
!   call Grad_Mod_For_Phi(grid, w % x, 3, var12y, .true.)  ! d2W/dxdz
!   call Grad_Mod_For_Phi(grid, w % y, 3, var12z, .true.)  ! d2W/dydz
!
!   do c = 1, grid % n_cells
!     Uxx = var1x(c)
!     Uyy = var1y(c)
!     Uzz = var1z(c)
!     Uxy = var2x(c)
!     Uxz = var2y(c)
!     Uyz = var2z(c)
!     Vxx = var9x(c)
!     Vyy = var9y(c)
!     Vzz = var9z(c)
!     Vxy = var10x(c)
!     Vxz = var10y(c)
!     Vyz = var10z(c)
!     Wxx = var11x(c)
!     Wyy = var11y(c)
!     Wzz = var11z(c)
!     Wxy = var12x(c)
!     Wxz = var12y(c)
!     Wyz = var12z(c)
!     dUdx= u % x(c) 
!     dUdy= u % y(c) 
!     dUdz= u % z(c) 
!     dVdx= v % x(c) 
!     dVdy= v % y(c) 
!     dVdz= v % z(c) 
!     dWdx= w % x(c) 
!     dWdy= w % y(c) 
!     dWdz= w % z(c) 
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
!     diss1(c) = duu_dx*Uxx + duv_dy*Uyy + duw_dz*Uzz   &
!              + Uxy*(duv_dx + duu_dy)                  &
!              + Uxz*(duw_dx + duu_dz)                  &
!              + Uyz*(duw_dy + duv_dz)                  &
!              + duv_dx*Vxx + dvv_dy*Vyy + dvw_dz*Vzz   &
!              + Vxy*(dvv_dx + duv_dy)                  &
!              + Vxz*(dvw_dx + duv_dz)                  &
!              + Vyz*(dvw_dy + dvv_dz)                  &
!              + duw_dx*Wxx + dvw_dy*Wyy + dww_dz*Wzz   &
!              + Wxy*(dvw_dx + duw_dy)                  &
!              + Wxz*(dww_dx + duw_dy)                  &
!              + Wyz*(dww_dy + dvw_dz)                  &
!              + 0.32 * kin%n(c) / eps%n(c)  *  &
!                (Uxx    * (duu_dx*dUdx  + duv_dx*dUdy  + duw_dx*dUdz ) + & 
!                 Uyy    * (duv_dy*dUdx  + dvv_dy*dUdy  + dvw_dy*dUdz ) + & 
!                 Uzz    * (duw_dz*dUdx  + dvw_dz*dUdy  + dww_dz*dUdz ) + & 
!                 Uxy    * (duu_dy*dUdx  + duv_dy*dUdy  + duw_dy*dUdz   + &
!                           duv_dx*dUdx  + dvv_dx*dUdy  + dvw_dx*dUdz ) + & 
!                 Uxz    * (duu_dz*dUdx  + duv_dz*dUdy  + duw_dz*dUdz   + &
!                           duw_dx*dUdx  + dvw_dx*dUdy  + dww_dx*dUdz ) + &
!                 Uyz    * (duv_dz*dUdx  + dvv_dz*dUdy  + dvw_dz*dUdz   + &
!                           duw_dy*dUdx  + dvw_dy*dUdy  + dww_dy*dUdz ) + &
!                 Vxx    * (duu_dx*dVdx  + duv_dx*dVdy  + duw_dx*dVdz ) + & 
!                 Vyy    * (duv_dy*dVdx  + dvv_dy*dVdy  + dvw_dy*dVdz ) + & 
!                 Vzz    * (duw_dz*dVdx  + dvw_dz*dVdy  + dww_dz*dVdz ) + & 
!                 Vxy    * (duu_dy*dVdx  + duv_dy*dVdy  + duw_dy*dVdz   + &
!                           duv_dx*dVdx  + dvv_dx*dVdy  + dvw_dx*dVdz ) + & 
!                 Vxz    * (duu_dz*dVdx  + duv_dz*dVdy  + duw_dz*dVdz   + &
!                           duw_dx*dVdx  + dvw_dx*dVdy  + dww_dx*dVdz ) + &
!                 Vyz    * (duv_dz*dVdx  + dvv_dz*dVdy  + dvw_dz*dVdz   + &
!                           duw_dy*dVdx  + dvw_dy*dVdy  + dww_dy*dVdz ) + &
!                 Wxx    * (duu_dx*dWdx  + duv_dx*dWdy  + duw_dx*dWdz ) + & 
!                 Wyy    * (duv_dy*dWdx  + dvv_dy*dWdy  + dvw_dy*dWdz ) + & 
!                 Wzz    * (duw_dz*dWdx  + dvw_dz*dWdy  + dww_dz*dWdz ) + & 
!                 Wxy    * (duu_dy*dWdx  + duv_dy*dWdy  + duw_dy*dWdz   + &
!                           duv_dx*dWdx  + dvv_dx*dWdy  + dvw_dx*dWdz ) + & 
!                 Wxz    * (duu_dz*dWdx  + duv_dz*dWdy  + duw_dz*dWdz   + &
!                           duw_dx*dWdx  + dvw_dx*dWdy  + dww_dx*dWdz ) + &
!                 Wyz    * (duv_dz*dWdx  + dvv_dz*dWdy  + dvw_dz*dWdz   + &
!                           duw_dy*dWdx  + dvw_dy*dWdy  + dww_dy*dWdz ))  
!     diss1(c) =  -2.0 * viscosity * diss1(c)
!   end do
! end if

  if(name_phi .eq. 'EPS') then
  do i=1,3
    if(i .eq. 1) then
      call Grad_Mod_For_Phi(grid, U % x, 1, ui_xx, .true.)  ! d2U/dxdx
      call Grad_Mod_For_Phi(grid, U % x, 2, ui_xy, .true.)  ! d2U/dxdy
      call Grad_Mod_For_Phi(grid, U % x, 3, ui_xz, .true.)  ! d2U/dxdz
      call Grad_Mod_For_Phi(grid, U % y, 2, ui_yy, .true.)  ! d2U/dydy
      call Grad_Mod_For_Phi(grid, U % y, 3, ui_yz, .true.)  ! d2U/dydz
      call Grad_Mod_For_Phi(grid, U % z, 3, ui_zz, .true.)  ! d2U/dzdz
    end if
    if(i .eq. 2) then
      call Grad_Mod_For_Phi(grid, V % x, 1, ui_xx, .true.)  ! d2V/dxdx
      call Grad_Mod_For_Phi(grid, V % x, 2, ui_xy, .true.)  ! d2V/dxdy
      call Grad_Mod_For_Phi(grid, V % x, 3, ui_xz, .true.)  ! d2V/dxdz
      call Grad_Mod_For_Phi(grid, V % y, 2, ui_yy, .true.)  ! d2V/dydy
      call Grad_Mod_For_Phi(grid, V % y, 3, ui_yz, .true.)  ! d2V/dydz
      call Grad_Mod_For_Phi(grid, V % z, 3, ui_zz, .true.)  ! d2V/dzdz
    end if
    if(i .eq. 3) then
      call Grad_Mod_For_Phi(grid, W % x, 1, ui_xx, .true.)  ! d2W/dxdx
      call Grad_Mod_For_Phi(grid, W % x, 2, ui_xy, .true.)  ! d2W/dxdy
      call Grad_Mod_For_Phi(grid, W % x, 3, ui_xz, .true.)  ! d2W/dxdz
      call Grad_Mod_For_Phi(grid, W % y, 2, ui_yy, .true.)  ! d2W/dydy
      call Grad_Mod_For_Phi(grid, W % y, 3, ui_yz, .true.)  ! d2W/dydz
      call Grad_Mod_For_Phi(grid, W % z, 3, ui_zz, .true.)  ! d2W/dzdz
    end if

    do c = 1, grid % n_cells
      if(i .eq. 1) then
        Uxx = ui_xx(c)
        Uxy = ui_xy(c)
        Uyx = Uxy
        Uxz = ui_xz(c)
        Uzx = Uxz
        Uyy = ui_yy(c)
        Uyz = ui_yz(c)
        Uzy = Uyz
        Uzz = ui_zz(c)
        diss1(c) =                                    &
                2.0*0.25*viscosity*kin%n(c)/eps_tot(c)  *  &
               (uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+  &
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+  &
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+  &
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+  &
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+  &
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+  &
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+  &
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+  &
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz))
      end if
      if(i .eq. 2) then
        Uxx = ui_xx(c)
        Uxy = ui_xy(c)
        Uyx = Uxy
        Uxz = ui_xz(c)
        Uzx = Uxz
        Uyy = ui_yy(c)
        Uyz = ui_yz(c)
        Uzy = Uyz
        Uzz = ui_zz(c)
        diss1(c) = diss1(c) +                         &
                2.0*0.25*viscosity*kin%n(c)/eps_tot(c)  *  &
               (uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+  &
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+  &
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+  &
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+  &
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+  &
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+  &
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+  &
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+  &
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz))
      end if
      if(i .eq. 3) then
        Uxx = ui_xx(c)
        Uxy = ui_xy(c)
        Uyx = Uxy
        Uxz = ui_xz(c)
        Uzx = Uxz
        Uyy = ui_yy(c)
        Uyz = ui_yz(c)
        Uzy = Uyz
        Uzz = ui_zz(c)
        diss1(c) = diss1(c) +                         &
                2.0*0.25*viscosity*kin%n(c)/eps_tot(c)  *  &
               (uu % n(c)*(Uxx*Uxx+Uxy*Uxy+Uxz*Uxz)+  &
                uv % n(c)*(Uxx*Uyx+Uxy*Uyy+Uxz*Uyz)+  &
                uw % n(c)*(Uxx*Uzx+Uxy*Uzy+Uxz*Uzz)+  &
                uv % n(c)*(Uyx*Uxx+Uyy*Uxy+Uyz*Uxz)+  &
                vv % n(c)*(Uyx*Uyx+Uyy*Uyy+Uyz*Uyz)+  &
                vw % n(c)*(Uyx*Uzx+Uyy*Uzy+Uyz*Uzz)+  &
                uw % n(c)*(Uzx*Uxx+Uzy*Uxy+Uzz*Uxz)+  &
                vw % n(c)*(Uzx*Uyx+Uzy*Uyy+Uzz*Uyz)+  &
                ww % n(c)*(Uzx*Uzx+Uzy*Uzy+Uzz*Uzz))
      end if
    end do
  end do  ! i
  end if                               

  call Grad_Mod_For_Phi(grid, l_scale, 1, l_sc_x,.true.) 
  call Grad_Mod_For_Phi(grid, l_scale, 2, l_sc_y,.true.) 
  call Grad_Mod_For_Phi(grid, l_scale, 3, l_sc_z,.true.) 

  r13 = ONE_THIRD
  r23 = TWO_THIRDS
  do  c = 1, grid % n_cells

    ! Epsilon over kinetic energy used almost 30 times in this loop
    eps_over_kin = eps % n(c) / kin%n(c)

    p_kin(c) = max( &
          - (  uu % n(c)*u % x(c) + uv % n(c)*u % y(c) + uw % n(c)*u % z(c)    &
             + uv % n(c)*v % x(c) + vv % n(c)*v % y(c) + vw % n(c)*v % z(c)    &
             + uw % n(c)*w % x(c) + vw % n(c)*w % y(c) + ww % n(c)*w % z(c)),  &
               1.0e-10)
  
    mag = max(0.0, sqrt(l_sc_x(c)**2 + l_sc_y(c)**2 + l_sc_z(c)**2), 1.0e-10)       

    n1 = l_sc_x(c) / mag 
    n2 = l_sc_y(c) / mag
    n3 = l_sc_z(c) / mag

    a11 = uu % n(c) / kin % n(c) - r23 
    a22 = vv % n(c) / kin % n(c) - r23
    a33 = ww % n(c) / kin % n(c) - r23
    a12 = uv % n(c) / kin % n(c)   
    a13 = uw % n(c) / kin % n(c)    
    a23 = vw % n(c) / kin % n(c)    
    
    s11 = u % x(c)
    s22 = v % y(c)
    s33 = w % z(c)
    s12 = 0.5*(u % y(c)+v % x(c))
    s13 = 0.5*(u % z(c)+w % x(c))
    s23 = 0.5*(v % z(c)+w % y(c))

    v11 = 0.0
    v22 = 0.0
    v33 = 0.0
    v12 = 0.5*(u % y(c)-v % x(c)) - omega_z
    v13 = 0.5*(u % z(c)-w % x(c)) + omega_y
    v23 = 0.5*(v % z(c)-w % y(c)) - omega_x

    aa2 = (a11**2)+(a22**2)+(a33**2)+2*((a12**2)+(a13**2)+(a23**2))

    aa3 = a11**3 + a22**3 + a33**3  &
        + 3*a12**2*(a11+a22)        &
        + 3*a13**2*(a11+a33)        &
        + 3*a23**2*(a22+a33)        &
        + 6*a12*a13*a23

    aa = 1.0 - (9.0/8.0)*(aa2-aa3)
    aa = max(aa,0.0)
    aa = min(aa,1.0)
 
    uu_nn = (uu % n(c)*n1*n1 + uv % n(c)*n1*n2 + uw % n(c)*n1*n3   &
           + uv % n(c)*n2*n1 + vv % n(c)*n2*n2 + vw % n(c)*n2*n3   &
           + uw % n(c)*n3*n1 + vw % n(c)*n3*n2 + ww % n(c)*n3*n3)

    a_mn_a_mn = a11*a11 + a22*a22 + a33*a33 + 2.0*(a12*a12+a13*a13+a23*a23)
    a_lk_s_lk = a11*s11 + a22*s22 + a33*s33 + 2.0*(a12*s12+a13*s13+a23*s23)
 
    ee=aa
    fss=1.0-(sqrt(aa) * ee**2)
    do icont=1,6
      eps11 = (1.0 - fss)*r23*eps %n(c) + fss * uu%n(c) * eps_over_kin     
      eps22 = (1.0 - fss)*r23*eps %n(c) + fss * vv%n(c) * eps_over_kin     
      eps33 = (1.0 - fss)*r23*eps %n(c) + fss * ww%n(c) * eps_over_kin     
      eps12 = fss * uv%n(c) * eps_over_kin     
      eps13 = fss * uw%n(c) * eps_over_kin     
      eps23 = fss * vw%n(c) * eps_over_kin      
      eps21 = eps12
      eps31 = eps13
      eps32 = eps23

      e11 = eps11 / eps%n(c) - r23
      e22 = eps22 / eps%n(c) - r23
      e33 = eps33 / eps%n(c) - r23
      e12 = eps12 / eps%n(c)
      e13 = eps13 / eps%n(c)
      e23 = eps23 / eps%n(c)
      e21 = e12
      e31 = e13
      e32 = e23
      E2  = (e11**2)+(e22**2)+(e33**2)+2*((e12**2)+(e13**2)+(e23**2))

      E3 =   e11**3 + e22**3 + e33**3  &
         + 3*e12**2*(e11+e22)          &
         + 3*e13**2*(e11+e33)          &
         + 3*e23**2*(e22+e33)          &
         + 6*e12*e13*e23

      ee = 1.0 - (9.0/8.0) * (E2 - E3)

      ee = max(ee,0.0)
      ee = min(ee,1.0)
      fss=1.0-(aa**0.5*ee**2.0)
    end do
     
    re_t  = (kin % n(c)**2)/(viscosity*eps_tot(c) + TINY)
    f_eps = 1.0 - ((c_2e-1.4)/c_2e)*exp(-(re_t/6.0)**2.0)
    ff2  = min((re_t/150)**1.5, 1.0)
    fd   = 1.0/(1.0+0.1*re_t)
    ff1  = min(0.6, aa2)
    CC   = 2.5*aa*ff1**0.25*ff2
    cc1  = CC+SQRT(aa)*(ee**2)
    cc2  = 0.8*SQRT(aa)
    c1w  = max((1.0-0.7*CC), 0.3)
    c2w  = min(aa,0.3)
    f_w  = min((kin%n(c)**1.5)/(2.5*eps_tot(c)*grid % wall_dist(c)),1.4)

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

    var2_11 = -cc2*(p11 - r23*p_kin(c))
    var2_22 = -cc2*(p22 - r23*p_kin(c))
    var2_33 = -cc2*(p33 - r23*p_kin(c))
    var2_12 = -cc2*p12
    var2_13 = -cc2*p13
    var2_23 = -cc2*p23

    phi2_nn =   var2_11*n1*n1  &
            + 2*var2_12*n1*n2  &
            + 2*var2_13*n1*n3  &
            +   var2_22*n2*n2  &
            + 2*var2_23*n2*n3  &
            +   var2_33*n3*n3  

    var1w_11 = c1w * f_w * eps_over_kin              &
             * (uu_nn-1.5*2.0*(   uu%n(c)*n1*n1*0.0  &
                                + uv%n(c)*n1*n2      &
                                + uw%n(c)*n1*n3))

    var1w_22 = c1w * f_w * eps_over_kin              &
             * (uu_nn-1.5*2.0*(   uv%n(c)*n2*n1      &
                                + vv%n(c)*n2*n2*0.0  &
                                + vw%n(c)*n2*n3))

    var1w_33 = c1w * f_w * eps_over_kin              &
             * (uu_nn-1.5*2.0*(   uw%n(c)*n3*n1      &
                                + vw%n(c)*n3*n2      &
                                + ww%n(c)*n3*n3*0.0))

    var1w_12 = c1w * f_w * eps_over_kin      &
             *(-1.5*(  uu%n(c)*n2*n1         &
                     + uv%n(c)*n2*n2*0.0     &
                     + uw%n(c)*n2*n3         &
                     + uv%n(c)*n1*n1*0.0     &
                     + vv%n(c)*n1*n2         &
                     + vw%n(c)*n1*n3)) 

    var1w_13 = c1w * f_w * eps_over_kin      &
             *(-1.5*(  uu%n(c)*n3*n1         &
                     + uv%n(c)*n3*n2         &
                     + uw%n(c)*n3*n3*0.0     &
                     + uw%n(c)*n1*n1*0.0     &
                     + vw%n(c)*n1*n2         &
                     + ww%n(c)*n1*n3))

    var1w_23 = c1w * f_w * eps_over_kin      &
             * (-1.5*(  uw%n(c)*n2*n1        &
                      + vw%n(c)*n2*n2*0.0    &
                      + ww%n(c)*n2*n3        &
                      + uv%n(c)*n3*n1        &
                      + vv%n(c)*n3*n2        &
                      + vw%n(c)*n3*n3)*0.0)

    var2w_11 = c2w * f_w * (phi2_nn - 1.5*2.0 * (  var2_11*n1*n1  &
                                                 + var2_12*n1*n2  &
                                                 + var2_13*n1*n3))

    var2w_22 = c2w * f_w * (phi2_nn - 1.5*2.0 * (  var2_12*n1*n2  &
                                                 + var2_22*n2*n2  &
                                                 + var2_23*n3*n2))

    var2w_33 = c2w * f_w * (phi2_nn - 1.5*2.0 * (  var2_13*n1*n3  &
                                                 + var2_23*n2*n3  &
                                                 + var2_33*n3*n3))

    var2w_12 = c2w*f_w*(-1.5*(var2_11*n2*n1 + var2_12*n2*n2 + var2_13*n2*n3 +  &
                              var2_12*n1*n1 + var2_22*n1*n2 + var2_23*n1*n3))
    var2w_13 = c2w*f_w*(-1.5*(var2_11*n3*n1 + var2_12*n3*n2 + var2_13*n3*n3 +  &
                              var2_13*n1*n1 + var2_23*n1*n2 + var2_33*n1*n3))
    var2w_23 = c2w*f_w*(-1.5*(var2_13*n2*n1 + var2_23*n2*n2 + var2_33*n2*n3 +  &
                              var2_12*n3*n1 + var2_22*n3*n2 + var2_23*n3*n3))

    !---------------!
    !   UU stress   !
    !---------------!
    if(name_phi .eq. 'UU') then

      b(c) = b(c) + (cc1*eps%n(c)*r23 + max(p11,      0.0)  &
                                      + max(var2_11,  0.0)  &
                                      + max(var1w_11, 0.0)  &
                                      + max(var2w_11, 0.0)) * grid % vol(c) 

      A % val(A % dia(c)) = A % val(A % dia(c))                         &
                          + (  cc1 *       eps_over_kin                 &
                             + c1w * f_w * eps_over_kin     *3.0*n1*n1  &
                             + fss *       eps_over_kin     ) * grid % vol(c) 

      A % val(A % dia(c)) = A % val(A % dia(c))          &
                          + (  max(-p11,      0.0)       &
                             + max(-var2_11,  0.0)       &
                             + max(-var1w_11, 0.0)       &
                             + max(-var2w_11, 0.0)       &
                             + (1.0-fss)*r23*eps%n(c) )  &
                          / max(uu%n(c),1.0e-10) * grid % vol(c) 

    !---------------!
    !   VV stress   !
    !---------------!
    else if(name_phi .eq. 'VV') then

      b(c) = b(c) + (cc1*eps%n(c)*r23 + max(p22,     0.0)  &
                                      + max(var2_22, 0.0)  &
                                      + max(var1w_22,0.0)  &
                                      + max(var2w_22,0.0)) * grid % vol(c) 

      A % val(A % dia(c)) = A % val(A % dia(c))                         &
                          + (  cc1 *       eps_over_kin                 &
                             + c1w * f_w * eps_over_kin     *3.0*n2*n2  &
                             + fss *       eps_over_kin     ) * grid % vol(c) 

      A % val(A % dia(c)) = A % val(A % dia(c))          &
                          + (  max(-p22,      0.0)       &
                             + max(-var2_22,  0.0)       &
                             + max(-var1w_22, 0.0)       &
                             + max(-var2w_22, 0.0)       &
                             + (1.0-fss)*r23*eps%n(c) )  &
                          / max(vv%n(c),1.0e-10) * grid % vol(c) 

    !---------------!
    !   WW stress   !
    !---------------!
    else if(name_phi .eq. 'WW') then

      b(c) = b(c) + (cc1*eps%n(c)*r23 + max(p33,     0.0)  &
                                      + max(var2_33, 0.0)  &
                                      + max(var1w_33,0.0)  &
                                      + max(var2w_33,0.0)) * grid % vol(c) 

      A % val(A % dia(c)) = A % val(A % dia(c))                         &
                          + (  cc1 *       eps_over_kin                 &
                             + c1w * f_w * eps_over_kin     *3.0*n3*n3  &
                             + fss *       eps_over_kin     ) * grid % vol(c) 

      A % val(A % dia(c)) = A % val(A % dia(c))          &
                          + (  max(-p33,      0.0)       &
                             + max(-var2_33,  0.0)       &
                             + max(-var1w_33, 0.0)       &
                             + max(-var2w_33, 0.0)       &
                             + (1.0-fss)*r23*eps%n(c) )  &
                          / max(ww%n(c),1.0e-10) * grid % vol(c) 

    !---------------!
    !   UV stress   !
    !---------------!
    else if(name_phi .eq. 'UV') then
      b(c) = b(c) + (p12 + var2_12 + var1w_12 + var2w_12) * grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))                                &
                          + (  cc1 *       eps_over_kin                        &
                             + c1w * f_w * eps_over_kin     *1.5*(n1*n1+n2*n2) &
                             + fss *       eps_over_kin      ) * grid % vol(c) 

    !---------------!
    !   UW stress   !
    !---------------!
    else if(name_phi .eq. 'UW') then
      b(c) = b(c) + (p13 + var2_13 + var1w_13 + var2w_13) * grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))                                &
                          + (  cc1 *       eps_over_kin                        &
                             + c1w * f_w * eps_over_kin     *1.5*(n1*n1+n3*n3) &
                             + fss *       eps_over_kin      ) * grid % vol(c) 

    !---------------!
    !   VW stress   !
    !---------------!
    else if(name_phi .eq. 'VW') then
      b(c) = b(c) + (p23 + var2_23 + var1w_23 + var2w_23) * grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))                                &
                          + (  cc1 *       eps_over_kin                        &
                             + c1w * f_w * eps_over_kin     *1.5*(n2*n2+n3*n3) &
                             + fss *       eps_over_kin      ) * grid % vol(c) 

    !----------------------!
    !   Epsilon equation   !
    !----------------------!
    else if(name_phi .eq. 'EPS') then 
      f_eps = 1.0 - ((c_2e-1.4)/c_2e) * exp(-(re_t/6.0)**2)
      eps_1 = 1.44 * p_kin(c) * eps_over_kin           
      eps_2 = c_2e * f_eps * eps_over_kin   
      b(c)  = b(c) + max(eps_1 + diss1(c), 0.0) * grid % vol(c) 
     
      A % val(A % dia(c)) =  A % val(A % dia(c)) + eps_2 * grid % vol(c)
    end if
  end do

  do c = 1, grid % n_cells
    kin_e(c) = sqrt( 0.5 * (uu % n(c) + vv % n(c) + ww % n(c)) )    
  end do 

  if(name_phi .eq. 'EPS') then
    call Grad_Mod_For_Phi(grid, kin_e, 1, kin_x, .true.)             ! dK/dx
    call Grad_Mod_For_Phi(grid, kin_e, 2, kin_y, .true.)             ! dK/dy
    call Grad_Mod_For_Phi(grid, kin_e, 3, kin_z, .true.)             ! dK/dz
    do c = 1, grid % n_cells
      re_t  = (kin % n(c)**2) / (viscosity*eps % n(c) + TINY)
      f_eps = 1.0 - ((c_2e-1.4)/c_2e) * exp(-(re_t/6.0)**2)
      b(c) = b(c) + (c_2e * f_eps * eps % n(c) / kin % n(c)                 &
                     * (viscosity*(kin_x(c)**2 + kin_y(c)**2 + kin_z(c)**2)))  &
                  * grid % vol(c)
    end do
  end if

  if(name_phi .eq. 'EPS') then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Calculate a values of dissipation  on wall
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER ) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          ! HOTFIXED IT (c->c2) - CHECK IT
          eps % n(c2) = viscosity*(kin_x(c2)**2 + kin_y(c2)**2 + kin_z(c2)**2)
        end if   ! end if of BC=wall
      end if    ! end if of c2<0
    end do
  end if

  end subroutine
