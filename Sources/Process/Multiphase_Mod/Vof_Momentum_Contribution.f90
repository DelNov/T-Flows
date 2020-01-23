  !============================================================================!
    subroutine Multiphase_Mod_Vof_Momentum_Contribution(mult, i, b)
  !----------------------------------------------------------------------------!
  !   Computes Surface tension and Gravity sources for Momentum Equation if    !
  !   a two-phase flow calculation is performed.                               !
  !----------------------------------------------------------------------------!
    implicit none
  !---------------------------------[Arguments]--------------------------------!
    type(Multiphase_Type), target :: mult
    real,                  target :: b(:)
    integer                       :: i
  !-----------------------------------[Locals]---------------------------------!
    type(Field_Type), pointer :: flow
    type(Grid_Type),  pointer :: grid
    type(Var_Type),   pointer :: vof
    integer                   :: s, c, c1, c2
    real                      :: dotprod, epsloc
    real,             pointer :: si(:), sj(:), sk(:)
    real                      :: corr_x, corr_y, corr_z
  !============================================================================!

    ! Take aliases
    flow   => mult % pnt_flow
    grid   => mult % pnt_grid
    vof    => mult % vof
    si     => grid % sx
    sj     => grid % sy
    sk     => grid % sz

    epsloc = epsilon(epsloc)

    ! Surface tension contribution
    if (mult % surface_tension > TINY) then

      select case(i)
        case(1)
          do c = 1, grid % n_cells
            b(c) = b(c) + mult % surface_tension * mult % curv(c)     &
                                                 * vof % x(c)         &
                                                 * grid % vol(c)
          end do
        case(2)
          do c = 1, grid % n_cells
            b(c) = b(c) + mult % surface_tension * mult % curv(c)     &
                                                 * vof % y(c)         &
                                                 * grid % vol(c)
          end do
        case(3)
          do c = 1, grid % n_cells
            b(c) = b(c) + mult % surface_tension * mult % curv(c)     &
                                                 * vof % z(c)         &
                                                 * grid % vol(c)
          end do

      end select

    end if

    ! Body force contribution (gravity)
    if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

      select case(i)

        case(1)

          if( .not. allocated(mult % body_fx)) then
            allocate(mult % body_fx(-grid % n_bnd_cells:grid % n_cells))
            allocate(mult % body_fy(-grid % n_bnd_cells:grid % n_cells))
            allocate(mult % body_fz(-grid % n_bnd_cells:grid % n_cells))
          end if

          mult % body_fx = 0.0; mult % body_fy = 0.0; mult % body_fz = 0.0

          do s = 1, grid % n_faces
            c1 = grid % faces_c(1,s)
            c2 = grid % faces_c(2,s)

            dotprod = ( (grid % xf(s) - grid % xc(c1)) * grav_x  &
                      + (grid % yf(s) - grid % yc(c1)) * grav_y  &
                      + (grid % zf(s) - grid % zc(c1)) * grav_z )

            mult % body_fx(c1) = mult % body_fx(c1)                      &
                               + flow % density_f(s) * si(s) * dotprod
            mult % body_fy(c1) = mult % body_fy(c1)                      &
                               + flow % density_f(s) * sj(s) * dotprod
            mult % body_fz(c1) = mult % body_fz(c1)                      &
                               + flow % density_f(s) * sk(s) * dotprod

            if (c2 > 0) then
              ! Correction for periodic faces:
              call Grid_Mod_Correction_Periodicity(grid, s,   &
                                                   corr_x, corr_y, corr_z)

              dotprod = ( (( grid % xf(s) + corr_x)      &
                           - grid % xc(c2)) * grav_x     &
                        + (( grid % yf(s) + corr_y)      &
                           - grid % yc(c2)) * grav_y     &
                        + (( grid % zf(s) + corr_z)      &
                           - grid % zc(c2)) * grav_z )

              mult % body_fx(c2) = mult % body_fx(c2)                      &
                                 - flow % density_f(s) * si(s) * dotprod
              mult % body_fy(c2) = mult % body_fy(c2)                      &
                                 - flow % density_f(s) * sj(s) * dotprod
              mult % body_fz(c2) = mult % body_fz(c2)                      &
                                 - flow % density_f(s) * sk(s) * dotprod

            end if
          end do

        call Grid_Mod_Exchange_Real(grid, mult % body_fx)
        call Grid_Mod_Exchange_Real(grid, mult % body_fy)
        call Grid_Mod_Exchange_Real(grid, mult % body_fz)

        do c = 1, grid % n_cells
          b(c) = b(c) + mult % body_fx(c)
        end do

      case(2)

        do c = 1, grid % n_cells
          b(c) = b(c) + mult % body_fy(c)
        end do

      case(3)

        do c = 1, grid % n_cells
          b(c) = b(c) + mult % body_fz(c)
        end do

    end select

  end if

  end subroutine
