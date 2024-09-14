!==============================================================================!
  subroutine Vis_T_Wale(Turb, Flow, Grid)
!------------------------------------------------------------------------------!
!  Compute SGS viscosity for 'LES' by using LES_WALE model.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  type(Field_Type), target :: Flow
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer                   :: c
  real                      :: s11, s22, s33,  s12, s13, s23,  s21, s31, s32
  real                      :: s11d,s22d,s33d, s12d,s13d,s23d, s21d,s31d,s32d
  real                      :: v11, v22, v33,  v12, v13, v23,  v21, v31, v32
  real                      :: sijd_sijd, shear2, vort2
  real, contiguous, pointer :: u_x(:), u_y(:), u_z(:)
  real, contiguous, pointer :: v_x(:), v_y(:), v_z(:)
  real, contiguous, pointer :: w_x(:), w_y(:), w_z(:)
!==============================================================================!

  call Work % Connect_Real_Cell(u_x, u_y, u_z,  &
                                v_x, v_y, v_z,  &
                                w_x, w_y, w_z)

  call Flow % Grad_Component(Grid, Flow % u % n, 1, u_x)
  call Flow % Grad_Component(Grid, Flow % u % n, 2, u_y)
  call Flow % Grad_Component(Grid, Flow % u % n, 3, u_z)

  call Flow % Grad_Component(Grid, Flow % v % n, 1, v_x)
  call Flow % Grad_Component(Grid, Flow % v % n, 2, v_y)
  call Flow % Grad_Component(Grid, Flow % v % n, 3, v_z)

  call Flow % Grad_Component(Grid, Flow % w % n, 1, w_x)
  call Flow % Grad_Component(Grid, Flow % w % n, 2, w_y)
  call Flow % Grad_Component(Grid, Flow % w % n, 3, w_z)

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   u_x,  &
  !$acc   v_y,  &
  !$acc   w_z,  &
  !$acc   v_x,  &
  !$acc   u_y,  &
  !$acc   u_z,  &
  !$acc   w_x,  &
  !$acc   w_y,  &
  !$acc   v_z,  &
  !$acc   flow_shear,  &
  !$acc   flow_vort,  &
  !$acc   turb_wale_v   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)
    s11 = u_x(c)
    s22 = v_y(c)
    s33 = w_z(c)
    s12 = 0.5*(v_x(c) + u_y(c))
    s13 = 0.5*(u_z(c) + w_x(c))
    s23 = 0.5*(w_y(c) + v_z(c))
    s21 = s12
    s31 = s13
    s32 = s23

    v11 = 0.0
    v22 = 0.0
    v33 = 0.0
    v12 = 0.5*(v_x(c) - u_y(c))
    v13 = 0.5*(u_z(c) - w_x(c))
    v23 = 0.5*(w_y(c) - v_z(c))
    v21 = -v12
    v31 = -v13
    v32 = -v23

    shear2 = 0.5 * flow_shear(c) * flow_shear(c)
    vort2  = 0.5 * flow_vort(c)  * flow_vort(c)

    s11d =  s11*s11 + s12*s12 + s13*s13   &
         - (v11*v11 + v12*v12 + v13*v13)  &
         - ONE_THIRD * (shear2 - vort2)

    s22d =  s12*s12 + s22*s22 + s23*s23   &
         - (v12*v12 + v22*v22 + v23*v23)  &
         - ONE_THIRD * (shear2 - vort2)

    s33d =  s13*s13 + s23*s23 + s33*s33   &
         - (v13*v13 + v23*v23 + v33*v33)  &
         - ONE_THIRD * (shear2 - vort2)

    s12d = s11*s12 + s12*s22 + s13*s32 + (v11*v12 + v12*v22 + v13*v32)
    s13d = s11*s13 + s12*s23 + s13*s33 + (v11*v13 + v12*v23 + v13*v33)
    s23d = s21*s13 + s22*s23 + s23*s33 + (v21*v13 + v22*v23 + v23*v33)

    s21d = s12d
    s31d = s13d
    s32d = s23d

    sijd_sijd = s11d*s11d + s22d*s22d + s33d*s33d  &
              + s12d*s12d + s13d*s13d + s23d*s23d

    ! Below: 1.5 == 3/2; 2.5 == 5/2;  1.25 == 5/4
    turb_wale_v(c) =  ( abs(sijd_sijd)**1.5 )    &
                   /  ( abs(shear2)**2.5 + abs(sijd_sijd)**1.25 + TINY)
  end do
  !$acc end parallel

  call Work % Disconnect_Real_Cell(u_x, u_y, u_z,  &
                                   v_x, v_y, v_z,  &
                                   w_x, w_y, w_z)

  end subroutine
