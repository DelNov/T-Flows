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
                      s_0        => r_cell_10
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
  integer                    :: irk, solv_iter, nod
  integer                    :: s, c, c1, c2, i_fac
  integer                    :: i_sub, n_sub, n_sub_param, n_iter
  integer                    :: t_steps_rk3
  real                       :: fs, norm_val, delta, surf_ng
  real                       :: res_grad, res_dir, sumwf
  real                       :: df, dcross, w_f
  real                       :: phi_1_p, phi_2_p, phi_f
  real                       :: part_1, part_2, part_3
  real                       :: d_tau
  real                       :: xmin, xmax, ymin, ymax, zmin, zmax
  real                       :: epsloc, eps_grid
!==============================================================================!

  ! Take aliases
  flow      => mult % pnt_flow
  grid      => flow % pnt_grid
  vof       => mult % vof
  dist_func => mult % dist_func
  a         => sol % a

  epsloc = epsilon(epsloc)

  if (mult % t_dist_scheme == RUNGE_KUTTA_3) then
    t_steps_rk3 = 3
  else
    t_steps_rk3 = 1
  end if

  d_tau = mult % c_tau * grid % min_vol ** ONE_THIRD

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

    s_0(c) = dist_func % n(c)                                         &
           / sqrt(  dist_func % n(c) ** 2                             &
                  + (norm2((/dist_func % x(c) ** 2,                   &
                             dist_func % y(c) ** 2,                   &
                             dist_func % z(c) ** 2/)) * delta) ** 2)
  end do

  solv_iter = 0
  n_iter = ceiling(mult % c_eps * grid % max_vol ** ONE_THIRD / d_tau)

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

      call Multiphase_Mod_Vof_Solver_Dist_Function_Cell_Loop(grid,        &
                                                             a,           &
                                                             dist_curr,   &
                                                             r_phi,       &
                                                             s_0,         &
                                                             grad_i,      &
                                                             grad_j,      &
                                                             grad_k)

      select case (irk)
      case(1)
        do c = 1, grid % n_cells
          dist_func % oo(c) = r_phi(c) * d_tau + dist_func % n(c)
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

    if (mult % t_dist_scheme .ne. RUNGE_KUTTA_3) then
      dist_func % n = dist_func % oo
    end if

  end do  ! end Solver loop

  call Grid_Mod_Exchange_Real(grid, dist_func % n)

  call Field_Mod_Grad_Variable(flow, dist_func)

  ! Find Heavyside function
  call Multiphase_Mod_Vof_Heavyside_Function(mult)

  end subroutine
