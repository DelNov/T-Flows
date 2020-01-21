!==============================================================================!
  subroutine Multiphase_Mod_Compute_Distance_Function(mult, sol, dt, n)
!------------------------------------------------------------------------------!
!   Calculates distance function using VOF as initialization. This is done     !
!   based on the paper of Dianat 2017, "A Coupled Level Set and Volume of      !
!   Fluid method for automotive exterior water management applications"        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: r_phi      => r_cell_05, &
                      grad_i     => r_cell_06, &
                      grad_j     => r_cell_07, &
                      grad_k     => r_cell_08, &
                      dist_curr  => r_cell_09, &
                      s_0        => r_cell_10, &
                      heavy_side => r_cell_11
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  real                          :: dt
  integer                       :: n    ! current temporal iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: vof
  type(Var_Type),    pointer :: dist_func
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  integer                    :: irk, solv_iter, nod
  integer                    :: s, c, c1, c2, i_fac
  integer                    :: i_sub, n_sub, n_sub_param, n_iter
  integer                    :: t_scheme, t_steps_rk3
  real                       :: fs, norm_val, delta, surf_ng
  real                       :: g_fx, g_fy, g_fz, eps_f
  real                       :: grad_x, grad_y, grad_z, grad_fn
  real                       :: dot_p1, dot_p2, res_grad, res_dir, sumwf
  real                       :: df, dcross, w_f
  real                       :: phi_1_p, phi_2_p, phi_f
  real                       :: part_1, part_2, part_3
  real                       :: d_tau, c_tau, c_eps, tol, r_fac(3)
  real                       :: xmin, xmax, ymin, ymax, zmin, zmax
  real                       :: corr_x, corr_y, corr_z
  real                       :: epsloc, eps_grid
  real, allocatable          :: a_i(:), b_i(:), c_i(:)
  character(len=80)          :: solver, name_t
!==============================================================================!

  ! Take aliases
  flow      => mult % pnt_flow
  grid      => flow % pnt_grid
  vof       => mult % vof
  dist_func => mult % dist_func

  a => sol % a
  b => sol % b % val

  epsloc = epsilon(epsloc)

  call Control_Mod_Factor_Fictitious_Time_Vof(c_tau)
  call Control_Mod_Factor_Number_Cells_Distance_Function_Vof(c_eps)
  call Control_Mod_Distance_Function_Time_Integration_Scheme(name_t)
  t_scheme = Numerics_Mod_Time_Integration_Scheme_Code(name_t)

  if (t_scheme == RK3) then
    t_steps_rk3 = 3
  else
    t_steps_rk3 = 1
  end if

  d_tau = c_tau * grid % min_vol ** ONE_THIRD

  !----------------------------------!
  !   Preliminar distance function   !
  !----------------------------------!
  do c = 1, grid % n_cells
    ! Find bounding box:
    xmin =  HUGE; ymin =  HUGE; zmin =  HUGE;
    xmax = -HUGE; ymax = -HUGE; zmax = -HUGE;

    do nod = 1, grid % cells_n_nodes(c)
      xmin = min(xmin, grid % xn(grid % cells_n(nod,c)))
      ymin = min(ymin, grid % yn(grid % cells_n(nod,c)))
      zmin = min(zmin, grid % zn(grid % cells_n(nod,c)))
      xmax = max(xmax, grid % xn(grid % cells_n(nod,c)))
      ymax = max(ymax, grid % yn(grid % cells_n(nod,c)))
      zmax = max(zmax, grid % zn(grid % cells_n(nod,c)))
    end do

    delta = 0.8 * maxval((/xmax - xmin, ymax - ymin, zmax - zmin/))

    dist_func % n(c) = (2.0 * vof % n(c) - 1.0) * delta
  end do

  !--------------------------!
  !   Values at boundaries   !
  !--------------------------!

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (c2 < 0) then
      dist_func % n(c2) = dist_func % n(c1)
    end if
  end do

  call Grid_Mod_Exchange_Real(grid, dist_func % n)

  call Field_Mod_Grad_Variable(flow, dist_func)

  !--------------------------------!
  !   Modified distance function   !
  !--------------------------------!

  do c = 1, grid % n_cells

    ! Find bounding box
    xmin =  HUGE; ymin =  HUGE; zmin =  HUGE;
    xmax = -HUGE; ymax = -HUGE; zmax = -HUGE;

    do nod = 1, grid % cells_n_nodes(c)
      xmin = min(xmin, grid % xn(grid % cells_n(nod,c)))
      ymin = min(ymin, grid % yn(grid % cells_n(nod,c)))
      zmin = min(zmin, grid % zn(grid % cells_n(nod,c)))
      xmax = max(xmax, grid % xn(grid % cells_n(nod,c)))
      ymax = max(ymax, grid % yn(grid % cells_n(nod,c)))
      zmax = max(zmax, grid % zn(grid % cells_n(nod,c)))
    end do

    delta = maxval((/xmax - xmin, ymax - ymin, zmax - zmin/))

    s_0(c) = dist_func % n(c)   &
           / sqrt(  dist_func % n(c) ** 2   &
                  + (norm2((/dist_func % x(c) ** 2,   &
                             dist_func % y(c) ** 2,   &
                             dist_func % z(c) ** 2/)) * delta) ** 2   &
                 )
  end do

  solv_iter = 0
  n_iter = ceiling(c_eps * grid % max_vol ** ONE_THIRD / d_tau)

  ! Solver loop
  do while (solv_iter < n_iter)

    ! Runge-Kutta loop
    do irk = 1, t_steps_rk3
      r_phi = 0.0

      select case (irk)
      case(1)
        call Field_Mod_Grad_Variable(flow, dist_func)
        grad_i = dist_func % x
        grad_j = dist_func % y
        grad_k = dist_func % z
        dist_curr = dist_func % n
      case(2)
        call Field_Mod_Grad_Component(flow, dist_func % oo, 1, grad_i)
        call Field_Mod_Grad_Component(flow, dist_func % oo, 2, grad_j)
        call Field_Mod_Grad_Component(flow, dist_func % oo, 3, grad_k)
        dist_curr = dist_func % oo
      case(3)
        call Field_Mod_Grad_Component(flow, dist_func % o, 1, grad_i)
        call Field_Mod_Grad_Component(flow, dist_func % o, 2, grad_j)
        call Field_Mod_Grad_Component(flow, dist_func % o, 3, grad_k)
        dist_curr = dist_func % o
      end select

      ! Loop on cells
      do c = 1, grid % n_cells

        ! Allocate a_i, b_i, c_i for faces
        if (allocated(a_i)) then
          deallocate(a_i, b_i, c_i)
        end if

        allocate( a_i(grid % cells_n_faces(c)),   &
                  b_i(grid % cells_n_faces(c)),   &
                  c_i(grid % cells_n_faces(c)) )

        a_i = 0.0; b_i = 0.0; c_i = 0.0

        ! loop on faces
        do i_fac = 1, grid % cells_n_faces(c)

          s = grid % cells_f(i_fac, c)
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          fs = grid % f(s)


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
            g_fx = fs * grad_i(c1) + (1.0 - fs) * grad_i(c2)
            g_fy = fs * grad_j(c1) + (1.0 - fs) * grad_j(c2)
            g_fz = fs * grad_k(c1) + (1.0 - fs) * grad_k(c2)

            dot_p1 = -dot_product((/g_fx, g_fy, g_fz/),                        &
                                 (/grid % sx(s), grid % sy(s), grid % sz(s)/))

            dot_p2 = -dot_product((/g_fx, g_fy, g_fz/),                        &
                                 (/grid % dx(s), grid % dy(s), grid % dz(s)/))


            grad_fn = dot_p1 + a % fc(s) * ( dist_curr(c1) - dist_curr(c2)   &
                                           - dot_p2 )

            ! Correction for periodicity
            corr_x = 0.0; corr_y = 0.0; corr_z = 0.0
            ! In x direction
            if( (abs(grid % per_x) > NANO) .and.                             &
                Math_Mod_Approx_Real(   abs(grid % dx(s))                    &
                                      + abs(grid % xc(c1) - grid % xc(c2)),  &
                                        grid % per_x) ) then
                corr_x = grid % per_x
            end if

            ! In y direction
            if( (abs(grid % per_y) > NANO) .and.                             &
                Math_Mod_Approx_Real(   abs(grid % dy(s))                    &
                                      + abs(grid % yc(c1) - grid % yc(c2)),  &
                                        grid % per_y) ) then
                corr_y = grid % per_y
            end if

            ! In z direction
            if( (abs(grid % per_z) > NANO) .and.                             &
                Math_Mod_Approx_Real(   abs(grid % dz(s))                    &
                                      + abs(grid % zc(c1) - grid % zc(c2)),  &
                                        grid % per_z) ) then
                corr_z = grid % per_z
            end if

            r_fac = -(/(grid % xf(s) + corr_x) - grid % xc(c),   &
                       (grid % yf(s) + corr_y) - grid % yc(c),   &
                       (grid % zf(s) + corr_z) - grid % zc(c)/)

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

          r_phi(c) = s_0(c) * ( 1.0 - sqrt( maxval(a_i)   &
                                          + maxval(b_i)   &
                                          + maxval(c_i) ) )

      end do  ! end loop cells

      select case (irk)
      case(1)
        do c = 1, grid % n_cells
          dist_func % oo(c) = r_phi(c) * d_tau   &
                            + dist_func % n(c)
        end do

        ! Update boundaries
        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if (c2 < 0) then
            dist_func % oo(c2) = dist_func % oo(c1)
          end if
        end do

      case(2)
        do c = 1, grid % n_cells
          dist_func % o(c) = 0.25 * r_phi(c) * d_tau   &
                           + 0.75 * dist_func % n(c)   &
                           + 0.25 * dist_func % oo(c)
        end do

        ! Update boundaries
        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if (c2 < 0) then
            dist_func % o(c2) = dist_func % o(c1)
          end if
        end do
      case(3)

        do c = 1, grid % n_cells
          dist_func % n(c) = TWO_THIRDS * r_phi(c) * d_tau   &
                           + ONE_THIRD * dist_func % n(c)    &
                           + TWO_THIRDS * dist_func % o(c)
        end do

        ! Update boundaries
        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if (c2 < 0) then
            dist_func % n(c2) = dist_func % n(c1)
          end if
        end do
      end select

    end do  ! end Runge-Kutta loop

    solv_iter = solv_iter + 1

    if (t_scheme .NE. RK3) then
      dist_func % n = dist_func % oo
    end if

    end do  ! end Solver loop

  call Grid_Mod_Exchange_Real(grid, dist_func % n)

  call Field_Mod_Grad_Variable(flow, dist_func)

  ! Find Heavyside function
  do c = 1, grid % n_cells

    eps_grid = c_eps * grid % vol(c) ** ONE_THIRD

    if (dist_func % n(c) > eps_grid) then
      heavy_side(c) = 1.0
    else if (dist_func % n(c) < -eps_grid) then
      heavy_side(c) = 0.0
    else
      heavy_side(c) = 0.5 * (1.0 + dist_func % n(c)/eps_grid        &
                           + 1.0 / PI * sin(PI * dist_func % n(c)   &
                           /eps_grid))
    end if

  end do

  ! Values at boundaries
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (c2 < 0) then
      heavy_side(c2) = heavy_side(c1)
    end if
  end do

  ! Find gradients
  call Field_Mod_Grad_Component(flow, heavy_side, 1, vof % x)
  call Field_Mod_Grad_Component(flow, heavy_side, 2, vof % y)
  call Field_Mod_Grad_Component(flow, heavy_side, 3, vof % z)

  end subroutine
