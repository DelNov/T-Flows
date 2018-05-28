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
  real    :: a11, a22, a33, a12, a13, a21, a31, a23, a32
  real    :: n1,n2,n3,AA2,AA3,AA,Ret,ff2,fd,FF1,CC,C1W,C2W,f_w,uu_nn
  real    :: e11,e12,e13,e21,e22,e23,e31,e32,e33
  real    :: Eps11,Eps12,Eps13,Eps21,Eps22,Eps23,Eps31,Eps32,Eps33
  real    :: Feps, phi2_nn
  real    :: fss,E2,E3,EE,CC1,CC2
  real    :: Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, Uyx, Uzx, Uzy
  real    :: r13, r23
  real    :: S11, S22, S33, S12, S13, S21, S31, S23, S32
  real    :: v11, v22, v33, v12, v13, v21, v31, v23, v32
  real    :: a_lk_s_lk, a_mn_a_mn
  real    :: VAR1w_11, VAR1w_22, VAR1w_33, VAR1w_12, VAR1w_13, VAR1w_23
  real    :: VAR2w_11, VAR2w_22, VAR2w_33, VAR2w_12, VAR2w_13, VAR2w_23
  real    :: VAR1_11, VAR1_22, VAR1_33, VAR1_12, VAR1_13, VAR1_23
  real    :: VAR2_11, VAR2_22, VAR2_33, VAR2_12, VAR2_13, VAR2_23
  real    :: P11, P22, P33, P12, P13, P23, Eps1, Eps2
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

  EE = 0.5
  AA = 0.5

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
    Eps_tot(c) = max(eps % n(c) + &
      0.5 * viscosity * (kin_xx(c) + kin_yy(c) + kin_zz(c)), TINY)
  end do

! !---------------------------------------------------!
! !   Below is one of versions of Hanjalic-Jakirlic   !
! !      model that required much more memory         !
! !---------------------------------------------------!
! if(name_phi .eq. "23") then
!   call Grad_Mod_For_Phi(grid, uu % n, 1, VAR3x, .true.) ! duu/dx  
!   call Grad_Mod_For_Phi(grid, uu % n, 2, VAR3y, .true.) ! duu/dy  
!   call Grad_Mod_For_Phi(grid, uu % n, 3, VAR3z, .true.) ! duu/dz  
!
!   call Grad_Mod_For_Phi(grid, vv % n, 1, VAR4x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, vv % n, 2, VAR4y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, vv % n, 3, VAR4z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, ww % n, 1, VAR5x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, ww % n, 2, VAR5y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, ww % n, 3, VAR5z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, uv % n, 1, VAR6x, .true.) ! duv/dx  
!   call Grad_Mod_For_Phi(grid, uv % n, 2, VAR6y, .true.) ! duv/dy  
!   call Grad_Mod_For_Phi(grid, uv % n, 3, VAR6z, .true.) ! duv/dz  
!
!   call Grad_Mod_For_Phi(grid, uw % n, 1, kin_x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, uw % n, 2, kin_y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, uw % n, 3, kin_z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, vw % n, 1, VAR8x, .true.) ! duw/dx  
!   call Grad_Mod_For_Phi(grid, vw % n, 2, VAR8y, .true.) ! duw/dy  
!   call Grad_Mod_For_Phi(grid, vw % n, 3, VAR8z, .true.) ! duw/dz  
!
!   call Grad_Mod_For_Phi(grid, u % x, 1, VAR1x, .true.)  ! d2U/dxdx
!   call Grad_Mod_For_Phi(grid, u % y, 2, VAR1y, .true.)  ! d2U/dydy
!   call Grad_Mod_For_Phi(grid, u % z, 3, VAR1z, .true.)  ! d2U/dzdz
!   call Grad_Mod_For_Phi(grid, u % x, 2, VAR2x, .true.)  ! d2U/dxdy
!   call Grad_Mod_For_Phi(grid, u % x, 3, VAR2y, .true.)  ! d2U/dxdz
!   call Grad_Mod_For_Phi(grid, u % y, 3, VAR2z, .true.)  ! d2U/dydz
!
!   call Grad_Mod_For_Phi(grid, v % x, 1, VAR9x, .true.)  ! d2V/dxdx
!   call Grad_Mod_For_Phi(grid, v % y, 2, VAR9y, .true.)  ! d2V/dydy
!   call Grad_Mod_For_Phi(grid, v % z, 3, VAR9z, .true.)  ! d2V/dzdz
!   call Grad_Mod_For_Phi(grid, v % x, 2, VAR10x, .true.)  ! d2V/dxdy
!   call Grad_Mod_For_Phi(grid, v % x, 3, VAR10y, .true.)  ! d2V/dxdz
!   call Grad_Mod_For_Phi(grid, v % y, 3, VAR10z, .true.)  ! d2V/dydz
!
!   call Grad_Mod_For_Phi(grid, w % x, 1, VAR11x, .true.)  ! d2W/dxdx
!   call Grad_Mod_For_Phi(grid, w % y, 2, VAR11y, .true.)  ! d2W/dydy
!   call Grad_Mod_For_Phi(grid, w % z, 3, VAR11z, .true.)  ! d2W/dzdz
!   call Grad_Mod_For_Phi(grid, w % x, 2, VAR12x, .true.)  ! d2W/dxdy
!   call Grad_Mod_For_Phi(grid, w % x, 3, VAR12y, .true.)  ! d2W/dxdz
!   call Grad_Mod_For_Phi(grid, w % y, 3, VAR12z, .true.)  ! d2W/dydz
!
!   do c = 1, grid % n_cells
!     Uxx = VAR1x(c)
!     Uyy = VAR1y(c)
!     Uzz = VAR1z(c)
!     Uxy = VAR2x(c)
!     Uxz = VAR2y(c)
!     Uyz = VAR2z(c)
!     Vxx = VAR9x(c)
!     Vyy = VAR9y(c)
!     Vzz = VAR9z(c)
!     Vxy = VAR10x(c)
!     Vxz = VAR10y(c)
!     Vyz = VAR10z(c)
!     Wxx = VAR11x(c)
!     Wyy = VAR11y(c)
!     Wzz = VAR11z(c)
!     Wxy = VAR12x(c)
!     Wxz = VAR12y(c)
!     Wyz = VAR12z(c)
!     dUdx= u % x(c) 
!     dUdy= u % y(c) 
!     dUdz= u % z(c) 
!     dVdx= v % x(c) 
!     dVdy= v % y(c) 
!     dVdz= v % z(c) 
!     dWdx= w % x(c) 
!     dWdy= w % y(c) 
!     dWdz= w % z(c) 
!     duu_dx = VAR3x(c)  
!     duu_dy = VAR3y(c)  
!     duu_dz = VAR3z(c)  
!     dvv_dx = VAR4x(c)  
!     dvv_dy = VAR4y(c)  
!     dvv_dz = VAR4z(c)  
!     dww_dx = VAR5x(c)  
!     dww_dy = VAR5y(c)  
!     dww_dz = VAR5z(c)  
!     duv_dx = VAR6x(c)  
!     duv_dy = VAR6y(c)  
!     duv_dz = VAR6z(c)  
!     duw_dx = kin_x(c)  
!     duw_dy = kin_y(c)  
!     duw_dz = kin_z(c)  
!     dvw_dx = VAR8x(c)  
!     dvw_dy = VAR8y(c)  
!     dvw_dz = VAR8z(c)  
!
!     Diss1(c) = duu_dx*Uxx + duv_dy*Uyy + duw_dz*Uzz   &
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
!     Diss1(c) =  -2.0 * viscosity * Diss1(c)
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
        Diss1(c) =                                    &
                2.0*0.25*viscosity*kin%n(c)/Eps_tot(c)  *  &
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
        Diss1(c) = Diss1(c) +                         &
                2.0*0.25*viscosity*kin%n(c)/Eps_tot(c)  *  &
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
        Diss1(c) = Diss1(c) +                         &
                2.0*0.25*viscosity*kin%n(c)/Eps_tot(c)  *  &
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
    p_kin(c) = max( &
          - (  uu % n(c)*u % x(c) + uv % n(c)*u % y(c) + uw % n(c)*u % z(c)    &
             + uv % n(c)*v % x(c) + vv % n(c)*v % y(c) + vw % n(c)*v % z(c)    &
             + uw % n(c)*w % x(c) + vw % n(c)*w % y(c) + ww % n(c)*w % z(c)),  &
               1.0e-10)
  
    mag = max(0.0, sqrt(l_sc_x(c)**2 + l_sc_y(c)**2 + l_sc_z(c)**2), 1.0e-10)       

    n1 = l_sc_x(c) / mag 
    n2 = l_sc_y(c) / mag
    n3 = l_sc_z(c) / mag

    a11 = uu % n(c)/kin % n(c) - r23 
    a22 = vv % n(c)/kin % n(c) - r23
    a33 = ww % n(c)/kin % n(c) - r23
    a12 = uv % n(c)/kin % n(c)   
    a21 = a12
    a13 = uw % n(c)/kin % n(c)    
    a31 = a13
    a23 = vw % n(c)/kin % n(c)    
    a32 = a23
    
    S11 = u % x(c)
    S22 = v % y(c)
    S33 = w % z(c)
    S12 = 0.5*(u % y(c)+v % x(c))
    S21 = S12
    S13 = 0.5*(u % z(c)+w % x(c))
    S31 = S13
    S23 = 0.5*(v % z(c)+w % y(c))
    S32 = S23

    v11 = 0.0
    v22 = 0.0
    v33 = 0.0
    v12 = 0.5*(u % y(c)-v % x(c)) - omega_z
    v21 = -v12 + omega_z
    v13 = 0.5*(u % z(c)-w % x(c)) + omega_y
    v31 = -v13 - omega_y
    v23 = 0.5*(v % z(c)-w % y(c)) - omega_x
    v32 = -v23 + omega_x

    AA2 = (a11**2)+(a22**2)+(a33**2)+2*((a12**2)+(a13**2)+(a23**2))

    AA3 = a11**3 + a22**3 + a33**3 +                 &
          3*a12**2*(a11+a22) + 3*a13**2*(a11+a33) +  &
          3*a23**2*(a22+a33) + 6*a12*a13*a23

    AA=1.0 - (9.0/8.0)*(AA2-AA3)
    AA=max(AA,0.0)
    AA=min(AA,1.0)
 
    uu_nn = (uu % n(c)*n1*n1+uv % n(c)*n1*n2+uw % n(c)*n1*n3   &
           + uv % n(c)*n2*n1+vv % n(c)*n2*n2+vw % n(c)*n2*n3   &
           + uw % n(c)*n3*n1+vw % n(c)*n3*n2+ww % n(c)*n3*n3)

    a_mn_a_mn = a11*a11 + a22*a22 + a33*a33 + 2.0*(a12*a12+a13*a13+a23*a23)
    a_lk_s_lk = a11*S11 + a22*S22 + a33*S33 + 2.0*(a12*S12+a13*S13+a23*S23)
 
    EE=AA
    fss=1.0-(sqrt(AA) * EE**2)
    do icont=1,6
      Eps11= (1.0 - fss)*r23*eps %n(c) + fss * uu%n(c)/kin%n(c) * eps%n(c)
      Eps22= (1.0 - fss)*r23*eps %n(c) + fss * vv%n(c)/kin%n(c) * eps%n(c)
      Eps33= (1.0 - fss)*r23*eps %n(c) + fss * ww%n(c)/kin%n(c) * eps%n(c)
      Eps12= fss * uv%n(c)/kin%n(c) * eps%n(c)
      Eps13= fss * uw%n(c)/kin%n(c) * eps%n(c)
      Eps23= fss * vw%n(c)/kin%n(c) * eps%n(c) 
      Eps21= Eps12
      Eps31= Eps13
      Eps32= Eps23

      e11= Eps11/eps%n(c) - r23
      e22= Eps22/eps%n(c) - r23
      e33= Eps33/eps%n(c) - r23
      e12= Eps12/eps%n(c)
      e13= Eps13/eps%n(c)
      e23= Eps23/eps%n(c)
      e21= e12
      e31= e13
      e32= e23
      E2=(e11**2)+(e22**2)+(e33**2)+2*((e12**2)+(e13**2)+(e23**2))

      E3= e11**3 + e22**3 + e33**3 + &
         3*e12**2*(e11+e22) + 3*e13**2*(e11+e33) +&
         3*e23**2*(e22+e33) + 6*e12*e13*e23

      EE=1.0 - (9.0/8.0)*(E2-E3)

      EE=max(EE,0.0)
      EE=min(EE,1.0)
      fss=1.0-(AA**0.5*EE**2.0)
    end do
     
    Ret= (kin % n(c)**2.)/(viscosity*eps_tot(c) + TINY)
    Feps = 1.0 - ((c_2e-1.4)/c_2e)*exp(-(Ret/6.0)**2.0)
    ff2 = min((Ret/150)**1.5, 1.0)
    fd=1.0/(1.0+0.1*Ret)
    FF1=min(0.6, AA2)
    CC=2.5*AA*FF1**0.25*ff2
    CC1=CC+SQRT(AA)*(EE**2)
    CC2=0.8*SQRT(AA)
    C1W=max((1.0-0.7*CC), 0.3)
    C2W=min(AA,0.3)
    f_w=min((kin%n(c)**1.5)/(2.5*Eps_tot(c)*grid % wall_dist(c)),1.4)

    P11 = - 2.0*(  uu % n(c) * u % x(c)     &
                 + uv % n(c) * u % y(c)     &
                 + uw % n(c) * u % z(c))    &
          - 2.0 * omega_y * 2.0 * uw % n(c)  &
          + 2.0 * omega_z * 2.0 * uv % n(c) 

    P22 = - 2.0*(  uv % n(c) * v % x(c)     &
                 + vv % n(c) * v % y(c)     &
                 + vw % n(c) * v % z(c))    &
          + 2.0 * omega_x * 2.0 * vw % n(c)  &
          - 2.0 * omega_z * 2.0 * uw % n(c) 

    P33 = - 2.0*(  uw % n(c) * w % x(c)     &
                 + vw % n(c) * w % y(c)     &
                 + ww % n(c) * w % z(c))    &
          - 2.0 * omega_x * 2.0 * vw % n(c)  &
          + 2.0 * omega_y * 2.0 * uw % n(c) 

    P12 = -(  uu % n(c) * v % x(c)      &
            + uv % n(c) * v % y(c)      &
            + uw % n(c) * v % z(c)      &
            + uv % n(c) * u % x(c)      &
            + vv % n(c) * u % y(c)      &
            + vw % n(c) * u % z(c))     &
            + 2.0 * omega_x * uw % n(c)  &
            - 2.0 * omega_y * vw % n(c)  &
            + 2.0 * omega_z * (vv % n(c) - uu % n(c)) 

    P13 = -(  uw % n(c)*u % x(c)                      &
            + vw % n(c)*u % y(c)                      &
            + ww % n(c)*u % z(c)                      &
            + uu % n(c)*w % x(c)                      &
            + uv % n(c)*w % y(c)                      &
            + uw % n(c)*w % z(c))                     &
            - 2.0 * omega_x * uv % n(c)                &
            - 2.0 * omega_y * (ww % n(c) - uu % n(c))  &
            + 2.0 * omega_z * vw % n(c) 

    P23 = -(  uv % n(c) * w % x(c)                    &
            + vv % n(c) * w % y(c)                    &
            + vw % n(c) * w % z(c)                    &
            + uw % n(c) * v % x(c)                    &
            + vw % n(c) * v % y(c)                    &
            + ww % n(c) * v % z(c))                   &
            - 2.0 * omega_x * (vw % n(c) - ww % n(c))  &
            + 2.0 * omega_y * uv % n(c)                &
            - 2.0 * omega_z * uw % n(c) 

    VAR1_11 = -CC1*eps%n(c)*a11 
    VAR1_22 = -CC1*eps%n(c)*a22 
    VAR1_33 = -CC1*eps%n(c)*a33 
    VAR1_12 = -CC1*eps%n(c)*a12 
    VAR1_13 = -CC1*eps%n(c)*a13 
    VAR1_23 = -CC1*eps%n(c)*a23 

    VAR2_11 = -CC2*(P11 - r23*p_kin(c))
    VAR2_22 = -CC2*(P22 - r23*p_kin(c))
    VAR2_33 = -CC2*(P33 - r23*p_kin(c))
    VAR2_12 = -CC2*P12
    VAR2_13 = -CC2*P13
    VAR2_23 = -CC2*P23

    phi2_nn = VAR2_11*n1*n1+2*VAR2_12*n1*n2+2*VAR2_13*n1*n3+VAR2_22*n2*n2+2*VAR2_23*n2*n3+VAR2_33*n3*n3  

    VAR1w_11 = C1W*f_w*eps%n(c)/kin%n(c)*(uu_nn-1.5*2.0*(uu%n(c)*n1*n1*0.0+uv%n(c)*n1*n2+uw%n(c)*n1*n3))
    VAR1w_22 = C1W*f_w*eps%n(c)/kin%n(c)*(uu_nn-1.5*2.0*(uv%n(c)*n2*n1+vv%n(c)*n2*n2*0.0+vw%n(c)*n2*n3))
    VAR1w_33 = C1W*f_w*eps%n(c)/kin%n(c)*(uu_nn-1.5*2.0*(uw%n(c)*n3*n1+vw%n(c)*n3*n2+ww%n(c)*n3*n3*0.0))
    VAR1w_12 = C1W*f_w*eps%n(c)/kin%n(c)*(-1.5*(uu%n(c)*n2*n1+uv%n(c)*n2*n2*0.0+uw%n(c)*n2*n3 +&
                                               uv%n(c)*n1*n1*0.0+vv%n(c)*n1*n2+vw%n(c)*n1*n3)) 
    VAR1w_13 = C1W*f_w*eps%n(c)/kin%n(c)*(-1.5*(uu%n(c)*n3*n1+uv%n(c)*n3*n2+uw%n(c)*n3*n3*0.0 +&
                                               uw%n(c)*n1*n1*0.0+vw%n(c)*n1*n2+ww%n(c)*n1*n3))
    VAR1w_23 = C1W*f_w*eps%n(c)/kin%n(c)*(-1.5*(uw%n(c)*n2*n1+vw%n(c)*n2*n2*0.0+ww%n(c)*n2*n3 +&
                                               uv%n(c)*n3*n1+vv%n(c)*n3*n2+vw%n(c)*n3*n3)*0.0)

    VAR2w_11 = C2W*f_w*(phi2_nn-1.5*2.0*(VAR2_11*n1*n1+VAR2_12*n1*n2+VAR2_13*n1*n3))
    VAR2w_22 = C2W*f_w*(phi2_nn-1.5*2.0*(VAR2_12*n1*n2+VAR2_22*n2*n2+VAR2_23*n3*n2))
    VAR2w_33 = C2W*f_w*(phi2_nn-1.5*2.0*(VAR2_13*n1*n3+VAR2_23*n2*n3+VAR2_33*n3*n3))
    VAR2w_12 = C2W*f_w*(-1.5*(VAR2_11*n2*n1+VAR2_12*n2*n2+VAR2_13*n2*n3 +&
                             VAR2_12*n1*n1+VAR2_22*n1*n2+VAR2_23*n1*n3))
    VAR2w_13 = C2W*f_w*(-1.5*(VAR2_11*n3*n1+VAR2_12*n3*n2+VAR2_13*n3*n3 +&
                             VAR2_13*n1*n1+VAR2_23*n1*n2+VAR2_33*n1*n3))
    VAR2w_23 = C2W*f_w*(-1.5*(VAR2_13*n2*n1+VAR2_23*n2*n2+VAR2_33*n2*n3 +&
                             VAR2_12*n3*n1+VAR2_22*n3*n2+VAR2_23*n3*n3))

    ! uu stress
    if(name_phi .eq. 'UU') then
!==============================================================================!
      b(c) = b(c) + (max(P11,0.0)+CC1*eps%n(c)*r23+max(VAR2_11,0.0)+max(VAR1w_11,0.0)+max(VAR2w_11,0.0))*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*eps%n(c)/kin%n(c)+C1W*f_w*eps%n(c)/kin%n(c)*3.0*n1*n1 + &
                      fss*eps%n(c)/kin%n(c))*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))+(max(-P11,0.0)+max(-VAR2_11,0.0)+max(-VAR1w_11,0.0)+max(-VAR2w_11,0.0) + &
                      (1.0-fss)*r23*eps%n(c))/max(uu%n(c),1.0e-10)*grid % vol(c) 
!==============================================================================!
    ! vv stress
    else if(name_phi .eq. 'VV') then
!==============================================================================!
      b(c) = b(c) + (max(P22,0.0)+CC1*eps%n(c)*r23+max(VAR2_22,0.0)+max(VAR1w_22,0.0)+max(VAR2w_22,0.0))*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*eps%n(c)/kin%n(c)+C1W*f_w*eps%n(c)/kin%n(c)*3.0*n2*n2 + &
                      fss*eps%n(c)/kin%n(c))*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))+(max(-P22,0.0)+max(-VAR2_22,0.0)+max(-VAR1w_22,0.0)+max(-VAR2w_22,0.0)+ &
                      (1.0-fss)*r23*eps%n(c))/max(vv%n(c),1.0e-10)*grid % vol(c) 
!==============================================================================!
    ! ww stress
    else if(name_phi .eq. 'WW') then
!==============================================================================!
      b(c) = b(c) + (max(P33,0.0)+CC1*eps%n(c)*r23+max(VAR2_33,0.0)+max(VAR1w_33,0.0)+max(VAR2w_33,0.0))*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*eps%n(c)/kin%n(c)+C1W*f_w*eps%n(c)/kin%n(c)*3.0*n3*n3 + &
                      fss*eps%n(c)/kin%n(c))*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c))+(max(-P33,0.0)+max(-VAR2_33,0.0)+max(-VAR1w_33,0.0)+max(-VAR2w_33,0.0)+ &
                      (1.0-fss)*r23*eps%n(c))/max(ww%n(c),1.0e-10)*grid % vol(c) 
!==============================================================================!
!==============================================================================!
    ! uv stress
    else if(name_phi .eq. 'UV') then
      b(c) = b(c) + (P12 + VAR2_12 + VAR1w_12 + VAR2w_12)*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*eps%n(c)/kin%n(c)+C1W*f_w*eps%n(c)/kin%n(c)*1.5*(n1*n1+n2*n2) + &
                      fss*eps%n(c)/kin%n(c))*grid % vol(c) 
!==============================================================================!
!==============================================================================!
    ! uw stress
    else if(name_phi .eq. 'UW') then
      b(c) = b(c) + (P13 + VAR2_13 + VAR1w_13 + VAR2w_13)*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*eps%n(c)/kin%n(c)+C1W*f_w*eps%n(c)/kin%n(c)*1.5*(n1*n1+n3*n3) + &
                      fss*eps%n(c)/kin%n(c))*grid % vol(c) 

!==============================================================================!
!==============================================================================!
    ! vw stress
    else if(name_phi .eq. 'VW') then
      b(c) = b(c) + (P23 + VAR2_23 + VAR1w_23 + VAR2w_23)*grid % vol(c) 
      A % val(A % dia(c)) = A % val(A % dia(c)) + (CC1*eps%n(c)/kin%n(c)+C1W*f_w*eps%n(c)/kin%n(c)*1.5*(n2*n2+n3*n3) + &
                      fss*eps%n(c)/kin%n(c))*grid % vol(c) 
!==============================================================================!
!==============================================================================!
    ! Epsilon equation
    else if(name_phi .eq. 'EPS') then 
      Feps = 1.0 - ((c_2e-1.4)/c_2e) * exp(-(Ret/6.0)**2)
      Eps1 = 1.44 * p_kin(c) * eps % n(c) / kin % n(c)
      Eps2 = c_2e*Feps*eps%n(c)/kin%n(c)
      b(c) = b(c) + max(Eps1 + Diss1(c),0.0)*grid % vol(c) 
     
      A % val(A % dia(c)) =  A % val(A % dia(c)) + Eps2*grid % vol(c)
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
      Ret  = (kin % n(c)**2) / (viscosity*eps % n(c) + TINY)
      Feps = 1.0 - ((c_2e-1.4)/c_2e) * exp(-(Ret/6.0)**2)
      b(c) = b(c) + (c_2e * Feps * eps % n(c) / kin % n(c)                 &
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