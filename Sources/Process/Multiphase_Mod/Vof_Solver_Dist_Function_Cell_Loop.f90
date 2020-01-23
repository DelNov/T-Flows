!==============================================================================!
  subroutine Multiphase_Mod_Vof_Solver_Dist_Function_Cell_Loop(grid,        &
                                                               a,           &
                                                               dist_curr,   &
                                                               r_phi,       &
                                                               s_0,         &
                                                               grad_i,      &
                                                               grad_j,      &
                                                               grad_k)
!------------------------------------------------------------------------------!
!    Loop on cells for calculation of the Distance function. This function     !
!    only reduces the size of Compute_Distance_Function
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)           :: grid
  type(Matrix_Type), target :: a
  real                      :: grad_i   (-grid % n_bnd_cells : grid % n_cells),&
                               grad_j   (-grid % n_bnd_cells : grid % n_cells),&
                               grad_k   (-grid % n_bnd_cells : grid % n_cells),&
                               dist_curr(-grid % n_bnd_cells : grid % n_cells),&
                               r_phi    (-grid % n_bnd_cells : grid % n_cells),&
                               s_0      (-grid % n_bnd_cells : grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer                   :: c, s, c1, c2
  integer                   :: i_fac
  real                      :: fs
  real                      :: dot_p1, dot_p2, grad_fn
  real                      :: g_fx, g_fy, g_fz
  real                      :: a_i(6), b_i(6), c_i(6), r_fac(3)
  real                      :: corr_x, corr_y, corr_z
!==============================================================================!

  do c = 1, grid % n_cells

    a_i = -HUGE; b_i = - HUGE; c_i = -HUGE

    ! Loop on faces
    do i_fac = 1, grid % cells_n_faces(c)

      s = grid % cells_f(i_fac, c)
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      corr_x = 0.0; corr_y = 0.0; corr_z = 0.0

      if (c1 .eq. c) then
        if (c2 > 0) then
          g_fx = fs * grad_i(c1) + (1.0 - fs) * grad_i(c2)
          g_fy = fs * grad_j(c1) + (1.0 - fs) * grad_j(c2)
          g_fz = fs * grad_k(c1) + (1.0 - fs) * grad_k(c2)
        else
          g_fx = grad_i(c1)
          g_fy = grad_j(c1)
          g_fz = grad_k(c1)
        end if

        dot_p1 = dot_product((/g_fx, g_fy, g_fz/),                         &
                             (/grid % sx(s), grid % sy(s), grid % sz(s)/))

        dot_p2 = dot_product((/g_fx, g_fy, g_fz/),                         &
                             (/grid % dx(s), grid % dy(s), grid % dz(s)/))

        grad_fn = dot_p1 + a % fc(s) * ( dist_curr(c2) - dist_curr(c1)   &
                                       - dot_p2 )

        r_fac = (/grid % xf(s) - grid % xc(c),   &
                  grid % yf(s) - grid % yc(c),   &
                  grid % zf(s) - grid % zc(c)/)
      else
        ! Correction for periodic faces:
        call Grid_Mod_Correction_Periodicity(grid, s,   &
                                             corr_x, corr_y, corr_z)

        g_fx = fs * grad_i(c1) + (1.0 - fs) * grad_i(c2)
        g_fy = fs * grad_j(c1) + (1.0 - fs) * grad_j(c2)
        g_fz = fs * grad_k(c1) + (1.0 - fs) * grad_k(c2)

        dot_p1 = -dot_product((/g_fx, g_fy, g_fz/),                        &
                             (/grid % sx(s), grid % sy(s), grid % sz(s)/))

        dot_p2 = -dot_product((/g_fx, g_fy, g_fz/),                        &
                             (/grid % dx(s), grid % dy(s), grid % dz(s)/))


        grad_fn = dot_p1 + a % fc(s) * ( dist_curr(c1) - dist_curr(c2)   &
                                       - dot_p2 )

        r_fac = -(/( grid % xf(s) + corr_x ) - grid % xc(c),   &
                   ( grid % yf(s) + corr_y ) - grid % yc(c),   &
                   ( grid % zf(s) + corr_z ) - grid % zc(c)/)

      end if

      grad_fn = grad_fn / grid % s(s)

      if( (dist_curr(c) > 0.0 .and. r_fac(1) > 0.0) .or.     &
          (dist_curr(c) < 0.0 .and. r_fac(1) < 0.0) ) then
        a_i(i_fac) = min(0.0, grad_fn * grid % sx(s) / grid % s(s)) ** 2
      else
        if( (dist_curr(c) > 0.0 .and. r_fac(1) < 0.0) .or.     &
            (dist_curr(c) < 0.0 .and. r_fac(1) > 0.0) ) then
          a_i(i_fac) = max(0.0, grad_fn * grid % sx(s) / grid % s(s)) ** 2
        end if
      end if

      if( (dist_curr(c) > 0.0 .and. r_fac(2) > 0.0) .or.     &
          (dist_curr(c) < 0.0 .and. r_fac(2) < 0.0) ) then
        b_i(i_fac) = min(0.0, grad_fn * grid % sy(s) / grid % s(s)) ** 2
      else
        if( (dist_curr(c) > 0.0 .and. r_fac(2) < 0.0) .or.     &
            (dist_curr(c) < 0.0 .and. r_fac(2) > 0.0) ) then
          b_i(i_fac) = max(0.0, grad_fn * grid % sy(s) / grid % s(s)) ** 2
        end if
      end if

      if( (dist_curr(c) > 0.0 .and. r_fac(3) > 0.0) .or.     &
          (dist_curr(c) < 0.0 .and. r_fac(3) < 0.0) ) then
        c_i(i_fac) = min(0.0, grad_fn * grid % sz(s) / grid % s(s)) ** 2
      else
        if( (dist_curr(c) > 0.0 .and. r_fac(3) < 0.0) .or.     &
            (dist_curr(c) < 0.0 .and. r_fac(3) > 0.0) ) then
          c_i(i_fac) = max(0.0, grad_fn * grid % sz(s) / grid % s(s)) ** 2
        end if
      end if

      end do  ! end loop faces

      r_phi(c) = s_0(c) * ( 1.0                                    &
               - sqrt( maxval(a_i(1:grid % cells_n_faces(c)))      &
                     + maxval(b_i(1:grid % cells_n_faces(c)))      &
                     + maxval(c_i(1:grid % cells_n_faces(c))) ) )

  end do  ! end loop cells

  end subroutine
