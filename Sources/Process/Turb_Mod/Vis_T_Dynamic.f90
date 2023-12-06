!==============================================================================!
  subroutine Vis_T_Dynamic(Turb)
!------------------------------------------------------------------------------!
!   Calculates Smagorinsky constant with dynamic procedure                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  integer                    :: c, c1, c2, s, sj, cj, nc, nb
  real                       :: u_a, v_a, w_a
  real                       :: uu_a, vv_a, ww_a, uv_a, uw_a, vw_a
  real                       :: m_11_a, m_22_a, m_33_a, m_12_a, m_13_a, m_23_a
  real                       :: l_11, l_22, l_33, l_12, l_13, l_23
  real                       :: m_11, m_22, m_33, m_12, m_13, m_23
  real                       :: m_dot_m, l_dot_m, l_g, l_f, vol_e
  real, contiguous,  pointer :: u_f(:), v_f(:), w_f(:)
  real, contiguous,  pointer :: uu_f(:), vv_f(:), ww_f(:)
  real, contiguous,  pointer :: uv_f(:), uw_f(:), vw_f(:)
  real, contiguous,  pointer :: m_11_f(:), m_22_f(:), m_33_f(:)
  real, contiguous,  pointer :: m_12_f(:), m_13_f(:), m_23_f(:)
  real, contiguous,  pointer :: shear_test(:)
!------------------------------------------------------------------------------!
!                                                                              !
!   C is derived from:    Lij_res = Lij_mod                                    !
!                                                                              !
!   res - resolved, mod - modeled                                              ! 
!                                                                              ! 
!   Lij_res = <UiUj> - <Ui><Uj>,                                               !
!                                                                              !
!   where <> denote test filter                                                !
!                                                                              !
!   Lij_mod = 2 * C * Mij, C = Csmag ** 2.0                                    !
!                                                                              !
!   where Mij is:  Mij = <delta**2.0>|<Sij>|<Sij> - <delta |Sij| Sij>          !
!                                                                              ! 
!   Finaly C is :                                                              ! 
!                                                                              !
!   C = 0.5 * Lij:Mij / Mij:Mij                                                !
!                                                                              !
!   aij : bij = a11 * b11 + a22 * b22 + a33 * b33                              !
!             + 2.0 * a12 * b12 + 2.0 a13 * b13 + 2.0 * a23 * b23              !
!                                                                              !
!==============================================================================!

  call Work % Connect_Real_Cell(u_f, v_f, w_f,                       &
                                uu_f, vv_f, ww_f, uv_f, uw_f, vw_f,  &
                                m_11_f, m_22_f, m_33_f,              &
                                m_12_f, m_13_f, m_23_f,              &
                                shear_test)

  ! Set all helping arrays to zero (important here)
  u_f   (:) = 0.0;  v_f   (:) = 0.0;  w_f   (:) = 0.0
  uu_f  (:) = 0.0;  vv_f  (:) = 0.0;  ww_f  (:) = 0.0
  uv_f  (:) = 0.0;  uw_f  (:) = 0.0;  vw_f  (:) = 0.0
  m_11_f(:) = 0.0;  m_22_f(:) = 0.0;  m_33_f(:) = 0.0
  m_12_f(:) = 0.0;  m_13_f(:) = 0.0;  m_23_f(:) = 0.0
  shear_test(:) = 0.0

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  nc   =  Grid % n_cells
  nb   =  Grid % n_bnd_cells
  call Flow % Alias_Momentum(u, v, w)

  call Grid % Exchange_Cells_Real(u % n)
  call Grid % Exchange_Cells_Real(v % n)
  call Grid % Exchange_Cells_Real(w % n)
  call Grid % Exchange_Cells_Real(Flow % shear)

  do c = Cells_In_Domain_And_Buffers()
    u_a   = 0.0
    v_a   = 0.0
    w_a   = 0.0
    vol_e = 0.0

    uu_a  = 0.0
    vv_a  = 0.0
    ww_a  = 0.0
    uv_a  = 0.0
    uw_a  = 0.0
    vw_a  = 0.0

    m_11_a = 0.0
    m_22_a = 0.0
    m_33_a = 0.0
    m_12_a = 0.0
    m_13_a = 0.0
    m_23_a = 0.0

    do sj = 1, Grid % cells_n_faces(c)  ! browse thrugh faces surrouind the cell
      s = Grid % cells_f(sj, c)         ! true face number
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      if(c2 .gt. 0) then
        if(c1 .eq. c) cj = c2
        if(c2 .eq. c) cj = c1

        ! Test velocities
        u_a = u_a + Grid % vol(cj) * u % n(cj)
        v_a = v_a + Grid % vol(cj) * v % n(cj)
        w_a = w_a + Grid % vol(cj) * w % n(cj)

        ! Test stresses
        uu_a = uu_a + Grid % vol(cj) * u % n(cj) * u % n(cj)
        vv_a = vv_a + Grid % vol(cj) * v % n(cj) * v % n(cj)
        ww_a = ww_a + Grid % vol(cj) * w % n(cj) * w % n(cj)
        uv_a = uv_a + Grid % vol(cj) * u % n(cj) * v % n(cj)
        uw_a = uw_a + Grid % vol(cj) * u % n(cj) * w % n(cj)
        vw_a = vw_a + Grid % vol(cj) * v % n(cj) * w % n(cj)

        ! Test Mija
        m_11_a = m_11_a + Grid % vol(cj) * Flow % shear(cj) * u % x(cj)
        m_22_a = m_22_a + Grid % vol(cj) * Flow % shear(cj) * v % y(cj)
        m_33_a = m_33_a + Grid % vol(cj) * Flow % shear(cj) * w % z(cj)
        m_12_a = m_12_a + Grid % vol(cj) * Flow % shear(cj)              &
               * 0.5 * ( u % y(cj) + v % x(cj) )
        m_13_a = m_13_a + Grid % vol(cj) * Flow % shear(cj)              &
               * 0.5 * ( u % z(cj) + w % x(cj) )
        m_23_a = m_23_a + Grid % vol(cj) * Flow % shear(cj)              &
               * 0.5 * ( v % z(cj) + w % y(cj) )

        ! Test volume 
        vol_e = vol_e + Grid % vol(cj) 
      end if
    end do

    ! Take into account influence of central cell within test molecule
    vol_e = vol_e + Grid % vol(c)

    u_a = u_a + Grid % vol(c) * u % n(c)
    v_a = v_a + Grid % vol(c) * v % n(c)
    w_a = w_a + Grid % vol(c) * w % n(c)

    uu_a = uu_a + Grid % vol(c) * u % n(c) * u % n(c)
    vv_a = vv_a + Grid % vol(c) * v % n(c) * v % n(c)
    ww_a = ww_a + Grid % vol(c) * w % n(c) * w % n(c)
    uv_a = uv_a + Grid % vol(c) * u % n(c) * v % n(c)
    uw_a = uw_a + Grid % vol(c) * u % n(c) * w % n(c)
    vw_a = vw_a + Grid % vol(c) * v % n(c) * w % n(c)

    m_11_a = m_11_a + Grid % vol(c) * Flow % shear(c) * u % x(c)
    m_22_a = m_22_a + Grid % vol(c) * Flow % shear(c) * v % y(c)
    m_33_a = m_33_a + Grid % vol(c) * Flow % shear(c) * w % z(c)
    m_12_a = m_12_a + Grid % vol(c) * Flow % shear(c) * .5*(u % y(c) + v % x(c))
    m_13_a = m_13_a + Grid % vol(c) * Flow % shear(c) * .5*(u % z(c) + w % x(c))
    m_23_a = m_23_a + Grid % vol(c) * Flow % shear(c) * .5*(v % z(c) + w % y(c))

    ! Now calculating test values
    u_f(c) = u_a / vol_e
    v_f(c) = v_a / vol_e
    w_f(c) = w_a / vol_e

    uu_f(c)  = uu_a / vol_e
    vv_f(c)  = vv_a / vol_e
    ww_f(c)  = ww_a / vol_e
    uv_f(c)  = uv_a / vol_e
    uw_f(c)  = uw_a / vol_e
    vw_f(c)  = vw_a / vol_e

    m_11_f(c) = m_11_a / vol_e
    m_22_f(c) = m_22_a / vol_e
    m_33_f(c) = m_33_a / vol_e
    m_12_f(c) = m_12_a / vol_e
    m_13_f(c) = m_13_a / vol_e
    m_23_f(c) = m_23_a / vol_e
  end do

  call Flow % Grad(Grid, u_f(-nb:nc), u % x,  &  ! dU/dx
                                      u % y,  &  ! dU/dy
                                      u % z)     ! dU/dz
  call Flow % Grad(Grid, v_f(-nb:nc), v % x,  &  ! dV/dx
                                      v % y,  &  ! dV/dy
                                      v % z)     ! dV/dz
  call Flow % Grad(Grid, w_f(-nb:nc), w % x,  &  ! dW/dx
                                      w % y,  &  ! dW/dy
                                      w % z)     ! dW/dz

  do c = Cells_In_Domain_And_Buffers()
    l_g  = Grid % vol(c)**ONE_THIRD
    l_f  = 2.0 * l_g

    shear_test(c) = sqrt(2.0*(  u % x(c)*u % x(c)                            &
                              + v % y(c)*v % y(c)                            &
                              + w % z(c)*w % z(c) +                          &
                         0.5*(v % z(c) + w % y(c))*(v % z(c) + w % y(c)) +   &
                         0.5*(u % z(c) + w % x(c))*(u % z(c) + w % x(c)) +   &
                         0.5*(v % x(c) + u % y(c))*(v % x(c) + u % y(c))))

    l_11 = uu_f(c) - u_f(c) * u_f(c)
    l_22 = vv_f(c) - v_f(c) * v_f(c)
    l_33 = ww_f(c) - w_f(c) * w_f(c)
    l_12 = uv_f(c) - u_f(c) * v_f(c)
    l_13 = uw_f(c) - u_f(c) * w_f(c)
    l_23 = vw_f(c) - v_f(c) * w_f(c)

    m_11 = l_f**2 * shear_test(c) * u % x(c) - l_g**2 * m_11_f(c)
    m_22 = l_f**2 * shear_test(c) * v % y(c) - l_g**2 * m_22_f(c)
    m_33 = l_f**2 * shear_test(c) * w % z(c) - l_g**2 * m_33_f(c)

    m_12 = l_f**2 * shear_test(c) * .5*(u % y(c)+v % x(c)) - l_g**2 * m_12_f(c)
    m_13 = l_f**2 * shear_test(c) * .5*(u % z(c)+w % x(c)) - l_g**2 * m_13_f(c)
    m_23 = l_f**2 * shear_test(c) * .5*(v % z(c)+w % y(c)) - l_g**2 * m_23_f(c)

    m_dot_m = m_11**2 + m_22**2 + m_33**2 + 2.0 * (m_12**2 + m_13**2 + m_23**2)

    l_dot_m =        l_11 * m_11 + l_22 * m_22 + l_33 * m_33   &
            + 2.0 * (l_12 * m_12 + l_13 * m_13 + l_23 * m_23)

    Turb % c_dyn(c)  =  -0.5 * l_dot_m / (m_dot_m + TINY)

    ! Set lower and upper limiter on c_dyn
    if(Turb % c_dyn(c) < 0.0) then
      Turb % c_dyn(c) = 0.0
    else if(Turb % c_dyn(c) > 0.04) then
      Turb % c_dyn(c) = 0.04
    end if
  end do

  call Work % Disconnect_Real_Cell(u_f, v_f, w_f,                       &
                                   uu_f, vv_f, ww_f, uv_f, uw_f, vw_f,  &
                                   m_11_f, m_22_f, m_33_f,              &
                                   m_12_f, m_13_f, m_23_f,              &
                                   shear_test)

  end subroutine
