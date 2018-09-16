!==============================================================================!
  subroutine Source_Ebm(grid, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for 'EBM'.                                                     !
!   Following paper "Recent progress in the development of                     !
!   the Elliptic Blending Reynolds-stress model"                               !
!   DOI: doi.org/10.1016/j.ijheatfluidflow.2014.09.002                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod

!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2
  real    :: prod_and_coriolis, phi_hom, phi_wall, Esor
  real    :: c_1e1
  real    :: stress
  real    :: eps_2_k, alpha3
  real    :: n1n1, n2n2, n3n3, n1n2, n1n3, n2n3, mag_f22
  real    :: b11, b22, b33, b12, b13, b23
  real    :: s11, s22, s33, s12, s13, s23
  real    :: v12, v13, v23
  real    :: b_kl_b_kl_sq, b_lm_s_lm, u_k_u_l_n_k_n_l, term_c3_1
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production P_ujj, p_kin  [m^2/s^3] | rate-of-strain  sij       [1/s]       !
!   dissipation   eps % n    [m^2/s^3] | turb. visc.     vis_t     [kg/(m*s)]  !
!   vel.-Pres. cor. phi      [m^2/s^3] | dyn visc.       viscosity [kg/(m*s)]  !
!   density       density    [kg/m^3]  | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol        [m^3]     | length          lf        [m]         !
!   left hand s.  A          [kg/s]    | kinematic viscosity       [m^2/s]     !
!   right hand s. for u_iu_j [m^2/s^3] | right hand s. for eps  b  [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
!   production : P_ij = - <u_i u_k> dU_j/dx_k - <u_j u_k> dU_j/dx_k            !
!   coriolis:    G_ij = - 2 omega_k ( eps_ikm <u_j u_m> - eps_jkm <u_i u_m> )  !
!------------------------------------------------------------------------------!

  call Time_And_Length_Scale(grid)

  do  c = 1, grid % n_cells

    ! no need to compute this for EPS -> can be improved
    kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
    ! P_k = 0.5 P_ii = - u_i u_k dU_i/dx_k
    p_kin(c) = max(-( uu % n(c) * u % x(c)  &
                     + uv % n(c) * u % y(c)  &
                     + uw % n(c) * u % z(c)  &
                     + uv % n(c) * v % x(c)  &
                     + vv % n(c) * v % y(c)  &
                     + vw % n(c) * v % z(c)  &
                     + uw % n(c) * w % x(c)  &
                     + vw % n(c) * w % y(c)  &
                     + ww % n(c) * w % z(c)), TINY)

    ! |df22/x_j|
    mag_f22 = max( f22 % x(c)**2. + f22 % y(c)**2. + f22 % z(c)**2., TINY )

    ! formula C.9 (n_i never appears individually, only as n_i * n_j)
    n1n1 = f22 % x(c)**2.        / mag_f22
    n2n2 = f22 % y(c)**2.        / mag_f22
    n3n3 = f22 % z(c)**2.        / mag_f22
    n1n2 = f22 % x(c)*f22 % y(c) / mag_f22
    n1n3 = f22 % x(c)*f22 % z(c) / mag_f22
    n2n3 = f22 % y(c)*f22 % z(c) / mag_f22

    ! frequently used expressions
    alpha3  = f22 % n(c)**3.
    eps_2_k = eps % n(c) / kin % n(c)

    ! formula C.4
    b11 = uu % n(c)/(2.*kin % n(c)) - ONE_THIRD 
    b22 = vv % n(c)/(2.*kin % n(c)) - ONE_THIRD
    b33 = ww % n(c)/(2.*kin % n(c)) - ONE_THIRD
    b12 = uv % n(c)/(2.*kin % n(c))   
    b13 = uw % n(c)/(2.*kin % n(c))    
    b23 = vw % n(c)/(2.*kin % n(c))

    ! formula C.5
    s11 = u % x(c) 
    s22 = v % y(c) 
    s33 = w % z(c) 
    s12 = 0.5*(u % y(c) + v % x(c))
    s13 = 0.5*(u % z(c) + w % x(c))
    s23 = 0.5*(v % z(c) + w % y(c))

    ! formula C.6
    v12 = 0.5*(u % y(c) - v % x(c)) - omega_z
    !v21 = -v12
    v13 = 0.5*(u % z(c) - w % x(c)) + omega_y
    !v31 = -v13
    v23 = 0.5*(v % z(c) - w % y(c)) - omega_x
    !v32 = -v23

    ! (C.3 1st term without "-" and *b_ij)
    term_c3_1 = g1*eps % n(c) + g1_star*p_kin(c)

    ! for formula C.3 (b_kl_b_kl never appears without sqrt)
    b_kl_b_kl_sq = sqrt(b11**2. + b22**2. + b33**2. &
      + 2.*(b12**2. + b13**2. + b23**2.))

    ! for formula C.3
    b_lm_s_lm = b11*s11 + b22*s22 + b33*s33 &
      + 2.*(b12*s12 + b13*s13 + b23*s23)

    ! for formula C.7
    u_k_u_l_n_k_n_l = uu % n(c)*n1n1 + vv % n(c)*n2n2 &
      + ww % n(c)*n3n3 + 2.*uv % n(c)*n1n2 + 2.*uw % n(c)*n1n3 &
      + 2.*vw % n(c)*n2n3

    !---------------!
    !   uu stress   !
    !---------------!
    if (name_phi .eq. 'UU') then
      ! limited stress
      stress = max(uu % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k * (                                 &
        2.*uu % n(c)*n1n1 + 2.*uv % n(c)*n1n2 + 2.*uw % n(c)*n1n3 &
        - 0.5*u_k_u_l_n_k_n_l*(n1n1 + 1.)                         &
        - stress) ! this extra term is substracted from A later
      
      ! formula C.3 (without C4 1st term)
      phi_hom = term_c3_1 * ONE_THIRD +                                  &
        ((g3 - g3_star*b_kl_b_kl_sq)*s11 +                               &
        g4*(2.*( b11*s11 + b12*s12 + b13*s13) - TWO_THIRDS*b_lm_s_lm ) + &
        g5*(2.*(           b12*v12 + b13*v13)))*kin % n(c)

      ! P_11 + G_11 (formula C.1)
      prod_and_coriolis = &
        -2.*(uu % n(c)*u % x(c) + uv % n(c)*u % y(c) + uw % n(c)*u % z(c)) &
        -2.*omega_y*2.*uw % n(c) + 2.*omega_z*2.*uv % n(c)

      ! left hand side (C.11 delta_ij)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * &
        TWO_THIRDS * alpha3 * eps % n(c) / stress

    !---------------!
    !   vv stress   !
    !---------------!
    elseif (name_phi .eq. 'VV') then

      ! limited stress
      stress = max(vv % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k * (                                 &
        2.*uv % n(c)*n1n2 + 2.*vv % n(c)*n2n2 + 2.*vw % n(c)*n2n3 &
        - 0.5*u_k_u_l_n_k_n_l*(n2n2 + 1.)                         &
        - stress) ! this extra term is substracted from A later

      ! formula C.3 (without C4 1st term)
      phi_hom = term_c3_1 * ONE_THIRD +                                  &
        ((g3 - g3_star*b_kl_b_kl_sq)*s22 +                               &
        g4*(2.*( b12*s12 + b22*s22 + b23*s23) - TWO_THIRDS*b_lm_s_lm ) + &
        g5*(2.*(-b12*v12           + b23*v23)))*kin % n(c)

      ! P_22 + G_22 (formula C.1)
      prod_and_coriolis = &
        -2.*(uv % n(c)*v % x(c) + vv % n(c)*v % y(c) + vw % n(c)*v % z(c)) &
        -2.*omega_x*2.*vw % n(c) + 2.*omega_z*2.*uw % n(c)

      ! left hand side (C.11 delta_ij)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * &
        TWO_THIRDS * alpha3 * eps % n(c) / stress
    !---------------!
    !   ww stress   !
    !---------------!
    elseif (name_phi .eq. 'WW') then

      ! limited stress
      stress = max(ww % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k * (                                 &
        2.*uw % n(c)*n1n3 + 2.*vw % n(c)*n2n3 + 2.*ww % n(c)*n3n3 &
        - 0.5*u_k_u_l_n_k_n_l*(n3n3 + 1.)                         &
        - stress) ! this extra term is substracted from A later

      ! formula C.3 (without C4 1st term)
      phi_hom = term_c3_1 * ONE_THIRD +                                 &
        ((g3 - g3_star*b_kl_b_kl_sq)*s33 +                              &
        g4*(2.*( b13*s13 + b23*s23 + b33*s33) - TWO_THIRDS*b_lm_s_lm) + &
        g5*(2.*(-b13*v13 - b23*v23           )))*kin % n(c)

      ! P_33 + G_33 (formula C.1)
      prod_and_coriolis = &
        -2.*(uw % n(c)*w % x(c) + vw % n(c)*w % y(c) + ww % n(c)*w % z(c)) &
        -2.*omega_x*2.*vw % n(c) + 2.*omega_y*2.*uw % n(c)

      ! left hand side (C.11 delta_ij)
      A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * &
        TWO_THIRDS * alpha3 * eps % n(c) / stress
    !---------------!
    !   uv stress   !
    !---------------!
    elseif (name_phi .eq. 'UV') then

      ! limited stress
      stress = max(uv % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k * (                           &
        uu % n(c)*n1n2 + uv % n(c)*n2n2 + uw % n(c)*n2n3 +  &
        uv % n(c)*n1n1 + vv % n(c)*n1n2 + vw % n(c)*n1n3    &
        - 0.5*u_k_u_l_n_k_n_l*n1n2                          &
        - stress) ! this extra term is substracted from A later

      ! formula C.3 (without C4 1st term)
      phi_hom = &
        ((g3 - g3_star*b_kl_b_kl_sq)*s12  +  &
        g4*( b11*s12 + b12*s22 + b13*s23  +  &
             b12*s11 + b22*s12 + b23*s13) +  &
        g5*(-b11*v12           + b13*v23  +  &
                       b22*v12 + b23*v13))*kin % n(c)

      ! P_12 + G_12 (formula C.1)
      prod_and_coriolis = &
        - uu % n(c)*v % x(c) - uv%n(c)*(v % y(c)+u % x(c)) - uw % n(c)*v % z(c)&
        - vv % n(c)*u % y(c) - vw % n(c)*u % z(c) &
        + 2.*omega_x*uw%n(c)-2.*omega_y*vw%n(c)+2.*omega_z*(vv%n(c)-uu%n(c))
    !---------------!
    !   uw stress   !
    !---------------!
    elseif (name_phi .eq. 'UW') then

      ! limited stress
      stress = max(uw % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k * (                               &
        uu % n(c)*n1n3 + uv % n(c)*n2n3 + uw % n(c)*n3n3+ &
        uw % n(c)*n1n1 + vw % n(c)*n1n2 + ww % n(c)*n1n3  &
        - 0.5*u_k_u_l_n_k_n_l*n1n3                           &
        - stress) ! this extra term is substracted from A later

      ! formula C.3 (without C4 1st term)
      phi_hom = &
        ((g3 - g3_star*b_kl_b_kl_sq)*s13  +  &
        g4*( b11*s13 + b12*s23 + b13*s33  +  &
             b13*s11 + b23*s12 + b33*s13) +  &
        g5*(-b11*v13 - b12*v23 +             &
                       b23*v12 + b33*v13))*kin % n(c)

      ! P_12 + G_12 (formula C.1)
      prod_and_coriolis = &
        - uu % n(c)*w % x(c) - uv%n(c)*w % y(c) - uw % n(c)*(w % z(c)+u % x(c))&
        - vw % n(c)*u % y(c) - ww % n(c)*u % z(c) &
        ! did not check
        - 2.*omega_x*uv%n(c) - 2.*omega_y*(ww%n(c)-uu%n(c)) + 2.*omega_z*vw%n(c)
    !---------------!
    !   vw stress   !
    !---------------!
    elseif (name_phi .eq. 'VW') then

      ! limited stress
      stress = max(vw % n(c), TINY)

      ! formula C.7
      phi_wall = - 5.*eps_2_k * (                          &
        uv % n(c)*n1n3 + vv % n(c)*n2n3 + vw % n(c)*n3n3 + &
        uw % n(c)*n1n2 + vw % n(c)*n2n2 + ww % n(c)*n2n3   &
        - 0.5*u_k_u_l_n_k_n_l*n2n3                         &
        - stress) ! this extra term is substracted from A later

      ! formula C.3 (without C4 1st term)
      phi_hom = &
        ((g3 - g3_star*b_kl_b_kl_sq)*s23  + &
        g4*( b12*s13 + b22*s23 + b23*s33  + &
             b13*s12 + b23*s22 + b33*s23) + &
        g5*(-b12*v13 - b22*v23              &
            -b13*v12 +           b33*v23))*kin % n(c)

      ! P_12 + G_12 (formula C.1)
      prod_and_coriolis = &
        - uu % n(c)*w % x(c) - uv%n(c)*w % y(c) - uw % n(c)*(w % z(c)+u % x(c))&
        - vw % n(c)*u % y(c) - ww % n(c)*u % z(c) &
        ! did not check
        - 2.*omega_x*(vv%n(c)-ww%n(c))+ 2.*omega_y*uv%n(c) - 2.*omega_z*uw%n(c)
    !-------------------------------------!
    !   repeating part for all stresses   !
    !-------------------------------------!
    ! formula C.1
    b(c) = b(c) + grid % vol(c) * ( & !
      max(prod_and_coriolis,0.) & ! P_ij + G_ij, if > 0.
      + (1. - alpha3)*phi_wall + alpha3*phi_hom & ! C.2
      )
    ! left hand side
    A % val(A % dia(c)) =  A % val(A % dia(c)) + grid % vol(c) * (  &
      + term_c3_1/(2.*kin % n(c)) & ! from C.3 and C.4 1st terms
      - min(prod_and_coriolis,0.)/stress +  & ! (P_ij + G_ij) / u_iu_j, if < 0
      6.*(1. - alpha3)*eps_2_k & ! substracting extra term from phi_wall
      )
    !---------!
    !   eps   !
    !---------!
    else if (name_phi .eq. 'EPS') then
      Esor = grid % vol(c)/max(t_scale(c),TINY)
      c_1e1 = c_1e * (1. + 0.1*(1.-alpha3)*p_kin(c)/(eps % n(c)+TINY))
      b(c) = b(c) + c_1e1*density*p_kin(c)*Esor

      ! Fill in a diagonal of coefficient matrix
      A % val(A % dia(c)) =  A % val(A % dia(c)) + c_2e*Esor*density
    end if 
  end do

  if (name_phi .eq. 'EPS') then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Calculate values of dissipation on wall
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          eps % n(c2) = (viscosity/density)*&
            (uu % n(c1) + vv % n(c1) + ww % n(c1))/grid % wall_dist(c1)**2.
        end if ! end if of BC=wall
      end if   ! end if of c2<0
    end do
  end if ! name_phi .eq. 'EPS'

  end subroutine
