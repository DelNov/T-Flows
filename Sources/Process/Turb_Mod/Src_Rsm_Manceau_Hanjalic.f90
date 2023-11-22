!==============================================================================!
  subroutine Src_Rsm_Manceau_Hanjalic(Turb, Sol, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for Re stresses and dissipation for                 !
!   RSM_MANCEAU_HANJALIC model                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
  character(len=*)          :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, zeta, f22
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c, s, c1, c2, nb, nc, reg
  real                       :: prod, diss, phi_hom, phi_wall, mag, phi_tot
  real                       :: b11, b22, b33, b12, b13, b21, b31, b23, b32
  real                       :: s11, s22, s33, s12, s13, s21, s31, s23, s32
  real                       :: v11, v22, v33, v12, v13, v21, v31, v23, v32
  real                       :: n1, n2, n3, b_mn_b_mn, b_lk_s_lk, ce_11
  real                       :: uu_nn, esor, diss_wall, diss_hom, kin_vis
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  nb   =  Grid % n_bnd_cells
  nc   =  Grid % n_cells
  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Sol % Alias_Native       (A, b)

  call Turb % Time_And_Length_Scale(Grid)

  call Flow % Grad_Variable(f22)

  do c = Cells_In_Domain()
    kin % n(c) = max(0.5*(  uu % n(c)  &
                          + vv % n(c)  &
                          + ww % n(c)), TINY)
    Turb % p_kin(c) = max(-(  uu % n(c)*u % x(c)  &
                            + uv % n(c)*u % y(c)  &
                            + uw % n(c)*u % z(c)  &
                            + uv % n(c)*v % x(c)  &
                            + vv % n(c)*v % y(c)  &
                            + vw % n(c)*v % z(c)  &
                            + uw % n(c)*w % x(c)  &
                            + vw % n(c)*w % y(c)  &
                            + ww % n(c)*w % z(c)), TINY)

    mag = max(TINY, sqrt(f22 % x(c)**2 + f22 % y(c)**2 + f22 % z(c)**2))
    n1 = f22 % x(c) / mag
    n2 = f22 % y(c) / mag
    n3 = f22 % z(c) / mag

    b11 = uu % n(c) / (2.0 * kin % n(c)) - ONE_THIRD
    b22 = vv % n(c) / (2.0 * kin % n(c)) - ONE_THIRD
    b33 = ww % n(c) / (2.0 * kin % n(c)) - ONE_THIRD
    b12 = uv % n(c) / (2.0 * kin % n(c))
    b21 = b12
    b13 = uw % n(c) / (2.0 * kin % n(c))
    b31 = b13
    b23 = vw % n(c) / (2.0 * kin % n(c))
    b32 = b23

    s11 = u % x(c) 
    s22 = v % y(c) 
    s33 = w % z(c) 
    s12 = 0.5*(u % y(c)+v % x(c)) 
    s21 = s12
    s13 = 0.5*(u % z(c)+w % x(c)) 
    s31 = s13 
    s23 = 0.5*(v % z(c)+w % y(c)) 
    s32 = s23 

    v11 = 0.0
    v22 = 0.0
    v33 = 0.0
    v12 =  0.5*(u % y(c)-v % x(c)) - Flow % omega_z
    v21 = -0.5*(u % y(c)-v % x(c)) + Flow % omega_z
    v13 =  0.5*(u % z(c)-w % x(c)) + Flow % omega_y
    v31 = -0.5*(u % z(c)-w % x(c)) - Flow % omega_y
    v23 =  0.5*(v % z(c)-w % y(c)) - Flow % omega_x
    v32 = -0.5*(v % z(c)-w % y(c)) + Flow % omega_x

    b_mn_b_mn = b11*b11 + b22*b22 + b33*b33 + 2.0*(b12*b12+b13*b13+b23*b23)
    b_lk_s_lk = b11*s11 + b22*s22 + b33*s33 + 2.0*(b12*s12+b13*s13+b23*s23)
    uu_nn     = (uu % n(c)*n1*n1+uv % n(c)*n1*n2+uw % n(c)*n1*n3 &
               + uv % n(c)*n2*n1+vv % n(c)*n2*n2+vw % n(c)*n2*n3 &
               + uw % n(c)*n3*n1+vw % n(c)*n3*n2+ww % n(c)*n3*n3)


    !---------------!
    !   UU stress   !
    !---------------!
    if(name_phi == 'UU') then
      phi_wall = -5.0 * eps % n(c) / kin % n(c) *      &
                 (-uu % n(c) + 2.0*uu % n(c)*n1*n1     &
                             + 2.0*uv % n(c)*n1*n2     &
                             + 2.0*uw % n(c)*n1*n3     &
                             - 0.5*(n1*n1+1.0)*uu_nn)

      phi_hom = - g1      * eps % n(c) * (-ONE_THIRD)                    &
                - g1_star * Turb % p_kin(c) * (-ONE_THIRD)               &
                + (g3-g3_star*sqrt(b_mn_b_mn)) * kin % n(c) * s11        &
                +  g4 * kin % n(c)*( 2.0*(b11*s11 + b12*s12 + b13*s13)   &
                                    -2.0/3.0 * b_lk_s_lk)                &
                +  g5 * kin % n(c)*( 2.0*(b11*v11 + b12*v12 + b13*v13))

      prod = - 2.0 * (  uu % n(c) * u % x(c)           &
                      + uv % n(c) * u % y(c)           &
                      + uw % n(c) * u % z(c))          &
             - 2.0 * Flow % omega_y * 2.0 * uw % n(c)  &
             + 2.0 * Flow % omega_z * 2.0 * uv % n(c)

      diss_wall = uu % n(c) / kin % n(c) * eps % n(c) 
      diss_hom  = 2.0/3.0 * eps % n(c)

      b(c) = b(c) + Flow % density(c) * (max(prod,0.0)               &
                  + (1.0-f22 % n(c)**2) * phi_wall                   &
                  +      f22 % n(c)**2  *(phi_hom)) * Grid % vol(c)

      A % val(A % dia(c)) = A % val(A % dia(c))                              &
        + Flow % density(c) * (max(-prod,0.0) / max(uu % n(c), TINY)         &
           + (1.0-f22 % n(c)**2) * 6.0 * eps % n(c)/kin % n(c)               &
           +      f22 % n(c)**2 *(diss_hom/max(uu % n(c), TINY)              &
                            + g1 * eps % n(c)           / (2.0*kin % n(c))   &
                            + g1_star * Turb % p_kin(c) / (2.0*kin % n(c)))  &
          ) * Grid % vol(c)

    !---------------!
    !   VV stress   !
    !---------------!
    else if(name_phi == 'VV') then
      phi_wall = -5.0 * eps % n(c) / kin % n(c) *      &
                 (-vv % n(c) + 2.0*uv % n(c)*n2*n1     &
                             + 2.0*vv % n(c)*n2*n2     &
                             + 2.0*vw % n(c)*n2*n3     &
                             - 0.5*(n2*n2+1.0)*uu_nn)

      phi_hom = - g1      * eps % n(c)      * (-ONE_THIRD)               &
                - g1_star * Turb % p_kin(c) * (-ONE_THIRD)               &
                + (g3-g3_star*sqrt(b_mn_b_mn)) * kin % n(c) * s22        &
                +  g4 * kin % n(c) *( 2.0*(b21*s21 + b22*s22 + b23*s23)  &
                                     -2.0/3.0 * b_lk_s_lk)               &
                +  g5 * kin % n(c) *( 2.0*(b21*v21 + b22*v22 + b23*v23))

      prod = - 2.0 * (  uv % n(c) * v % x(c)           &
                      + vv % n(c) * v % y(c)           &
                      + vw % n(c) * v % z(c))          &
             + 2.0 * Flow % omega_x * 2.0 *vw % n(c)   &
             - 2.0 * Flow % omega_z * 2.0 *uw % n(c)


      phi_tot = (1.0-f22 % n(c)**2) * phi_wall  &
              +      f22 % n(c)**2  * phi_hom

      diss_wall = vv % n(c)/kin % n(c) * eps % n(c) 
      diss_hom  = 2.0/3.0 * eps % n(c)

      b(c) = b(c) + Flow % density(c) * (max(prod,0.0)               &
                  + (1.0-f22 % n(c)**2) * phi_wall                   &
                  +      f22 % n(c)**2  *(phi_hom)) * Grid % vol(c)

      A % val(A % dia(c)) = A % val(A % dia(c))                              &
        + Flow % density(c) * (max(-prod,0.0) / max(vv % n(c), TINY)         &
           + (1.0-f22 % n(c)**2) * 6.0 * eps % n(c)/kin % n(c)               &
           +      f22 % n(c)**2 *(diss_hom/max(vv % n(c), TINY)              &
                            + g1 * eps % n(c)           / (2.0*kin % n(c))   &
                            + g1_star * Turb % p_kin(c) / (2.0*kin % n(c)))  &
          ) * Grid % vol(c)

    !---------------!
    !   WW stress   !
    !---------------!
    else if(name_phi == 'WW') then
      phi_wall = -5.0 * eps % n(c) / kin % n(c) *      &
                 (-ww % n(c) + 2.0*uw % n(c)*n3*n1     &
                             + 2.0*vw % n(c)*n3*n2     &
                             + 2.0*ww % n(c)*n3*n3     &
                             - 0.5*(n3*n3+1.0)*uu_nn)

      phi_hom = - g1      * eps % n(c)      * (-ONE_THIRD)               &
                - g1_star * Turb % p_kin(c) * (-ONE_THIRD)               &
                + (g3-g3_star*sqrt(b_mn_b_mn)) * kin % n(c) * s33        &
                +  g4 * kin % n(c) *( 2.0*(b31*s31 + b32*s32 + b33*s33)  &
                                     -2.0/3.0 * b_lk_s_lk)               &
                +  g5 * kin % n(c) *( 2.0*(b31*v31 + b32*v32 + b33*v33))

      prod = - 2.0 * (  uw % n(c) * w % x(c)           &
                      + vw % n(c) * w % y(c)           &
                      + ww % n(c) * w % z(c))          &
             - 2.0 * Flow % omega_x * 2.0 * vw % n(c)  &
             + 2.0 * Flow % omega_y * 2.0 * uw % n(c)

      phi_tot = (1.0-f22 % n(c)**2) * phi_wall  &
              +      f22 % n(c)**2  * phi_hom

      diss_wall = ww % n(c)/kin % n(c) * eps % n(c) 
      diss_hom  = 2.0/3.0 * eps % n(c)

      b(c) = b(c) + Flow % density(c) * (max(prod,0.0)               &
                  + (1.0-f22 % n(c)**2) * phi_wall                   &
                  +      f22 % n(c)**2  *(phi_hom)) * Grid % vol(c)

      A % val(A % dia(c)) = A % val(A % dia(c))                              &
        + Flow % density(c) * (max(-prod,0.0) / max(ww % n(c), TINY)         &
           + (1.0-f22 % n(c)**2) * 6.0 * eps % n(c)/kin % n(c)               &
           +      f22 % n(c)**2 *(diss_hom/max(ww % n(c), TINY)              &
                            + g1 * eps % n(c)           / (2.0*kin % n(c))   &
                            + g1_star * Turb % p_kin(c) / (2.0*kin % n(c)))  &
          ) * Grid % vol(c)

    !---------------!
    !   UV stress   !
    !---------------!
    else if(name_phi == 'UV') then
      phi_wall = -5.0 * eps % n(c) / kin % n(c) *  &
                 (-uv % n(c) + uu % n(c)*n2*n1     &
                             + uv % n(c)*n2*n2     &
                             + uw % n(c)*n2*n3     &
                             + uv % n(c)*n1*n1     &
                             + vv % n(c)*n1*n2     &
                             + vw % n(c)*n1*n3     &
                             - 0.5*(n1*n2)*uu_nn)

      phi_hom = (g3-g3_star*sqrt(b_mn_b_mn)) * kin % n(c) * s12  &
              +  g4*kin % n(c)*(  b11*s21 + b12*s22 + b13*s23    &
                                + b21*s11 + b22*s12 + b23*s13)   &
              +  g5*kin % n(c)*(  b11*v21 + b12*v22 + b13*v23    &
                                + b21*v11 + b22*v12 + b23*v13)

      prod = -(  uu % n(c) * v % x(c)                           &
               + uw % n(c) * v % z(c)                           &
               + uv % n(c) *(v % y(c) + u % x(c))               &
               + vv % n(c) * u % y(c)                           &
               + vw % n(c) * u % z(c))                          &
               + 2.0 * Flow % omega_x * uw % n(c)               &
               - 2.0 * Flow % omega_y * vw % n(c)               &
               + 2.0 * Flow % omega_z* (vv % n(c) - uu % n(c))

      phi_tot = (1.0-f22 % n(c)**2) * phi_wall  &
              +      f22 % n(c)**2  * phi_hom

      diss_wall = uv % n(c)/kin % n(c) * eps % n(c) 

      b(c) = b(c) + Flow % density(c) * (prod                        &
                  + (1.0-f22 % n(c)**2) * phi_wall                   &
                  +      f22 % n(c)**2  *(phi_hom)) * Grid % vol(c)

      A % val(A % dia(c)) =  A % val(A % dia(c))                              &
        + Flow % density(c) * ((1.0 - f22 % n(c)**2) * 6.0 * eps % n(c)       &
                            / kin % n(c)                                      &
                + f22 % n(c)**2 *(  g1 * eps % n(c) / (2.0 * kin % n(c))      &
                        + g1_star * Turb % p_kin(c) / (2.0 * kin % n(c)))     &
          ) * Grid % vol(c)

    !---------------!
    !   UW stress   !
    !---------------!
    else if(name_phi == 'UW') then
      phi_wall = -5.0 * eps % n(c) / kin % n(c) *  &
                 (-uw % n(c) + uu % n(c)*n3*n1     &
                             + uv % n(c)*n3*n2     &
                             + uw % n(c)*n3*n3     &
                             + uw % n(c)*n1*n1     &
                             + vw % n(c)*n1*n2     &
                             + ww % n(c)*n1*n3     &
                             - 0.5*(n1*n3)*uu_nn)

      phi_hom = (g3-g3_star*sqrt(b_mn_b_mn)) * kin % n(c) * s13  &
              +  g4*kin % n(c)*(  b11*s31 + b12*s32 + b13*s33    &
                                + b31*s11 + b32*s12 + b33*s13)   &
              +  g5*kin % n(c)*(  b11*v31 + b12*v32 + b13*v33    &
                                + b31*v11 + b32*v12 + b33*v13)

      prod = -(  uu % n(c) * w % x(c)                           &
               + uv % n(c) * w % y(c)                           &
               + uw % n(c) *(w % z(c) + u % x(c))               &
               + vw % n(c) * u % y(c)                           &
               + ww % n(c) * u % z(c))                          &
               - 2.0 * Flow % omega_x * uv % n(c)               &
               - 2.0 * Flow % omega_y *(ww % n(c) - uu % n(c))  &
               + 2.0 * Flow % omega_z * vw % n(c)

      phi_tot = (1.0-f22 % n(c)**2) * phi_wall  &
                   + f22 % n(c)**2  * phi_hom

      diss_wall = uw % n(c)/kin % n(c) * eps % n(c) 

      b(c) = b(c) + Flow % density(c) * (prod                        &
                  + (1.0-f22 % n(c)**2) * phi_wall                   &
                  +      f22 % n(c)**2  *(phi_hom)) * Grid % vol(c)

      A % val(A % dia(c)) =  A % val(A % dia(c))                               &
        + Flow % density(c) * ((1.0 - f22 % n(c)**2) * 6.0 * eps % n(c)    &
                            / kin % n(c)                                   &
                + f22 % n(c)**2  *(  g1 * eps % n(c) / (2.0 * kin % n(c))      &
                         + g1_star * Turb % p_kin(c) / (2.0 * kin % n(c)))     &
          ) * Grid % vol(c)

    !---------------!
    !   VW stress   !
    !---------------!
    else if(name_phi == 'VW') then
      phi_wall = -5.0 * eps % n(c) / kin % n(c) *  &
                 (-vw % n(c) + uv % n(c)*n3*n1     &
                             + vv % n(c)*n3*n2     &
                             + vw % n(c)*n3*n3     &
                             + uw % n(c)*n2*n1     &
                             + vw % n(c)*n2*n2     &
                             + ww % n(c)*n2*n3     &
                             - 0.5*(n2*n3)*uu_nn)

      phi_hom = (g3-g3_star*sqrt(b_mn_b_mn)) * kin % n(c) * s23  &
              +  g4*kin % n(c)*(  b21*s31 + b22*s32 + b23*s33    &
                                + b31*s21 + b32*s22 + b33*s23)   &
              +  g5*kin % n(c)*(  b21*v31 + b22*v32 + b23*v33    &
                                + b31*v21 + b32*v22 + b33*v23)

      prod = -(  uv % n(c) * w % x(c)                           &
               + vv % n(c) * w % y(c)                           &
               + vw % n(c) *(w % z(c) + v % y(c))               &
               + uw % n(c) * v % x(c)                           &
               + ww % n(c) * v % z(c))                          &
               - 2.0 * Flow % omega_x *(vw % n(c) - ww % n(c))  &
               + 2.0 * Flow % omega_y * uv % n(c)               &
               - 2.0 * Flow % omega_z * uw % n(c)

      phi_tot = (1.0-f22 % n(c)**2) * phi_wall  &
              +      f22 % n(c)**2  * phi_hom

      diss = (1.0 - f22 % n(c)**2) * vw % n(c) / kin % n(c) * eps % n(c)

      b(c) = b(c) + Flow % density(c) * (prod                        &
                  + (1.0-f22 % n(c)**2) * phi_wall                   &
                  +      f22 % n(c)**2  *(phi_hom)) * Grid % vol(c)

      A % val(A % dia(c)) =  A % val(A % dia(c))                               &
        + Flow % density(c) * ((1.0 - f22 % n(c)**2) * 6.0 * eps % n(c)  &
                           / kin % n(c)                                  &
        + f22 % n(c)**2  *(  g1 * eps % n(c) / (2.0 * kin % n(c))              &
                 + g1_star * Turb % p_kin(c) / (2.0 * kin % n(c)))             &
                     ) * Grid % vol(c)

    !----------------------!
    !   Epsilon equation   !
    !----------------------!
    else if(name_phi == 'EPS') then
      esor = Grid % vol(c) / max(Turb % t_scale(c), TINY)

      ce_11 = c_1e * (1.0 + 0.065*(1.0 - f22 % n(c)**3)  &
                   * Turb % p_kin(c) / eps % n(c))
      b(c) = b(c) + Flow % density(c) * ce_11 * Turb % p_kin(c) * esor

      ! Fill in a diagonal of coefficient matrix
      A % val(A % dia(c)) =  A % val(A % dia(c))  &
                          + c_2e * esor * Flow % density(c)
    end if
  end do
  call Grid % Exchange_Cells_Real(kin % n(-nb:nc))
  call Grid % Exchange_Cells_Real(Turb % p_kin(-nb:nc))

  if(name_phi == 'EPS') then
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. WALL .or.  &
         Grid % region % type(reg) .eq. WALLFL) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          ! Calculate a values of dissipation on wall
          kin_vis = Flow % viscosity(c1) / Flow % density(c1)
          eps % n(c2) = kin_vis*(uu % n(c1) + vv % n(c1) + ww % n(c1))&
                      / Grid % wall_dist(c1)**2
        end do  ! faces
      end if    ! boundary condition type
    end do      ! regions
  end if

  end subroutine
