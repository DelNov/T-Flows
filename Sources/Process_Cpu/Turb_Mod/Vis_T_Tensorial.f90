!==============================================================================!
  subroutine Vis_T_Tensorial(Turb)
!------------------------------------------------------------------------------!
!   Calculates the non-linear, tensorial "viscosity", and the associated SGS   !
!   stress tensor.                                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c
  real                      :: ixp, iyp, izp, ixyp, ixzp, iyzp
  real                      :: du_dx, du_dy, du_dz
  real                      :: dv_dx, dv_dy, dv_dz
  real                      :: dw_dx, dw_dy, dw_dz
  real                      :: sgv_x , sgv_y , sgv_z
  real                      :: sgv_xy, sgv_yx, sgv_xz
  real                      :: sgv_zx, sgv_yz, sgv_zy
!------------------------------------------------------------------------------!
!                                                                              !
!   nii_ki = (1/2*Vol) * I'_kh * dUi/dxh                                       !
!                                                                              !
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  ! Initialize the tensors
  Turb % ten_turb_11 = 0.0
  Turb % ten_turb_12 = 0.0
  Turb % ten_turb_13 = 0.0
  Turb % ten_turb_21 = 0.0
  Turb % ten_turb_22 = 0.0
  Turb % ten_turb_23 = 0.0
  Turb % ten_turb_31 = 0.0
  Turb % ten_turb_32 = 0.0
  Turb % ten_turb_33 = 0.0

  Turb % tau_11 = 0.0
  Turb % tau_12 = 0.0
  Turb % tau_13 = 0.0
  Turb % tau_21 = 0.0
  Turb % tau_22 = 0.0
  Turb % tau_23 = 0.0
  Turb % tau_31 = 0.0
  Turb % tau_32 = 0.0
  Turb % tau_33 = 0.0

  do c = Cells_In_Domain_And_Buffers()
    ixp = ( 0.5 * (  Grid % ixx (c)  &
                   - Grid % iyy (c)  &
                   - Grid % izz (c) ) )     ! I'xx         (I'11)
    iyp = ( 0.5 * (- Grid % ixx (c)  &
                   + Grid % iyy (c)  &
                   - Grid % izz (c) ) )     ! I'yy         (I'22)
    izp = ( 0.5 * (- Grid % ixx (c)  &
                   - Grid % iyy (c)  &
                   + Grid % izz (c) ) )     ! I'zz         (I'33)
    ixyp = Grid % ixy (c)              ! I'xy = I'yx  (I'12 = I'21)
    ixzp = Grid % ixz (c)              ! I'xz = I'zx  (I'13 = I'31)
    iyzp = Grid % iyz (c)              ! I'yz = I'zy  (I'23 = I'32)

    du_dx = u % x(c)  ! dU/dx
    du_dy = u % y(c)  ! dU/dy
    du_dz = u % z(c)  ! dU/dz
    dv_dx = v % x(c)  ! dV/dx
    dv_dy = v % y(c)  ! dV/dy
    dv_dz = v % z(c)  ! dV/dz
    dw_dx = w % x(c)  ! dW/dx
    dw_dy = w % y(c)  ! dW/dy
    dw_dz = w % z(c)  ! dW/dz

    ! Construct the non-linear "viscosity" tensor
    Turb % ten_turb_11 (c) = (0.5 / Grid % vol(c)) * ( ixp  * du_dx  &
                                                     + ixyp * du_dy  &
                                                     + ixzp * du_dz )
    Turb % ten_turb_22 (c) = (0.5 / Grid % vol(c)) * ( ixyp * dv_dx  &
                                                     + iyp  * dv_dy  &
                                                     + iyzp * dv_dz )
    Turb % ten_turb_33 (c) = (0.5 / Grid % vol(c)) * ( ixzp * dw_dx  &
                                                     + iyzp * dw_dy  &
                                                     + izp  * dw_dz )
    Turb % ten_turb_12 (c) = (0.5 / Grid % vol(c)) * ( ixp  * dv_dx  &
                                                     + ixyp * dv_dy  &
                                                     + ixzp * dv_dz )
    Turb % ten_turb_21 (c) = (0.5 / Grid % vol(c)) * ( ixyp * du_dx  &
                                                     + iyp  * du_dy  &
                                                     + iyzp * du_dz )
    Turb % ten_turb_13 (c) = (0.5 / Grid % vol(c)) * ( ixp  * dw_dx  &
                                                     + ixyp * dw_dy  &
                                                     + ixzp * dw_dz )
    Turb % ten_turb_31 (c) = (0.5 / Grid % vol(c)) * ( ixzp * du_dx  &
                                                     + iyzp * du_dy  &
                                                     + izp  * du_dz )
    Turb % ten_turb_23 (c) = (0.5 / Grid % vol(c)) * ( ixyp * dw_dx  &
                                                     + iyp  * dw_dy  &
                                                     + iyzp * dw_dz )
    Turb % ten_turb_32 (c) = (0.5 / Grid % vol(c)) * ( ixzp * dv_dx  &
                                                     + iyzp * dv_dy  &
                                                     + izp  * dv_dz )

    ! Since we're at it, we may as well construct the subgrid stress tensor too
    sgv_x  = Turb % ten_turb_11 (c)
    sgv_y  = Turb % ten_turb_22 (c)
    sgv_z  = Turb % ten_turb_33 (c)
    sgv_xy = Turb % ten_turb_12 (c)
    sgv_yx = Turb % ten_turb_21 (c)
    sgv_xz = Turb % ten_turb_13 (c)
    sgv_zx = Turb % ten_turb_31 (c)
    sgv_yz = Turb % ten_turb_23 (c)
    sgv_zy = Turb % ten_turb_32 (c)

    ! Now actually build the thing
    Turb % tau_11 (c)    = (- sgv_x  * du_dx - sgv_x  * du_dx)  &
                         + (- sgv_yx * du_dy - sgv_yx * du_dy)  &
                         + (- sgv_zx * du_dz - sgv_zx * du_dz)

    Turb % tau_22 (c)    = (- sgv_xy * dv_dx - sgv_xy * dv_dx)  &
                         + (- sgv_y  * dv_dy - sgv_y  * dv_dy)  &
                         + (- sgv_zy * dv_dz - sgv_zy * dv_dz)

    Turb % tau_33 (c)    = (- sgv_xz * dw_dx - sgv_xz * dw_dx)  &
                         + (- sgv_yz * dw_dy - sgv_yz * dw_dy)  &
                         + (- sgv_z  * dw_dz - sgv_z  * dw_dz)

    Turb % tau_12 (c)    = (- sgv_xy * du_dx - sgv_x  * dv_dx)  &
                         + (- sgv_y  * du_dy - sgv_yx * dv_dy)  &
                         + (- sgv_zy * du_dz - sgv_zx * dv_dz)

    Turb % tau_21 (c)    = (- sgv_x  * dv_dx - sgv_xy * du_dx)  &
                         + (- sgv_yx * dv_dy - sgv_y  * du_dy)  &
                         + (- sgv_zx * dv_dz - sgv_zy * du_dz)

    Turb % tau_13 (c)    = (- sgv_xz * du_dx - sgv_x  * dw_dx)  &
                         + (- sgv_yz * du_dy - sgv_yx * dw_dy)  &
                         + (- sgv_z  * du_dz - sgv_zx * dw_dz)

    Turb % tau_31 (c)    = (- sgv_x  * dw_dx - sgv_xz * du_dx)  &
                         + (- sgv_yx * dw_dy - sgv_yz * du_dy)  &
                         + (- sgv_zx * dw_dz - sgv_z  * du_dz)

    Turb % tau_23 (c)    = (- sgv_xz * dv_dx - sgv_xy * dw_dx)  &
                         + (- sgv_yz * dv_dy - sgv_y  * dw_dy)  &
                         + (- sgv_z  * dv_dz - sgv_zy * dw_dz)

    Turb % tau_32 (c)    = (- sgv_xy * dw_dx - sgv_xz * dv_dx)  &
                         + (- sgv_y  * dw_dy - sgv_yz * dv_dy)  &
                         + (- sgv_zy * dw_dz - sgv_z  * dv_dz)

  end do

  ! Call Save_Debug_Vtu to check if the tensors are calculated properly
  ! print *, 'SGS non-linear viscosity and stress tensors calculated'
  ! call Grid % Save_Debug_Vtu('new-nii-sgs-ij',                  &
  !                        tensor_cell = (/Turb % ten_turb_11,    &
  !                                        Turb % ten_turb_12,    &
  !                                        Turb % ten_turb_13,    &
  !                                        Turb % ten_turb_21,    &
  !                                        Turb % ten_turb_22,    &
  !                                        Turb % ten_turb_23,    &
  !                                        Turb % ten_turb_31,    &
  !                                        Turb % ten_turb_32,    &
  !                                        Turb % ten_turb_33/),  &
  !                        tensor_comp = 9,                       &
  !                        tensor_name = 'Nonlinear Turbulent Viscosity Tensor')

  ! call Grid % Save_Debug_Vtu('new-tau-sgs-ij',             &
  !                        tensor_cell = (/Turb % tau_11,    &
  !                                        Turb % tau_12,    &
  !                                        Turb % tau_13,    &
  !                                        Turb % tau_21,    &
  !                                        Turb % tau_22,    &
  !                                        Turb % tau_23,    &
  !                                        Turb % tau_31,    &
  !                                        Turb % tau_32,    &
  !                                        Turb % tau_33/),  &
  !                        tensor_comp = 9,                  &
  !                        tensor_name = 'Nonlinear Stress Tensor')

  end subroutine
