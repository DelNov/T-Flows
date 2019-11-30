  !============================================================================!
    subroutine Multiphase_Mod_Vof_Momentum_Contribution(mult, i, b)
  !----------------------------------------------------------------------------!
  !   Correct fluxes on pressure equation                                      !
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
  !============================================================================!

    ! Take aliases
    flow   => mult % pnt_flow
    grid   => mult % pnt_grid
    vof    => mult % vof
    si     => grid % sx
    sj     => grid % sy
    sk     => grid % sz

    epsloc = epsilon(epsloc)

    if (surface_tension > TINY) then

      select case(i)
        case(1)
          do c = 1, grid % n_cells
            b(c) = b(c) + surface_tension * mult % curv(c)     &
                                          * mult % vof % x(c)  &
                                          * grid % vol(c)
          end do
        case(2)
          do c = 1, grid % n_cells
            b(c) = b(c) + surface_tension * mult % curv(c)     &
                                          * mult % vof % y(c)  &
                                          * grid % vol(c)
          end do
        case(3)
          do c = 1, grid % n_cells
            b(c) = b(c) + surface_tension * mult % curv(c)     &
                                          * mult % vof % z(c)  &
                                          * grid % vol(c)
          end do

      end select

    end if


    if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

      select case(i)

        case(1)

          if(.not.allocated(body_fx)) then
            allocate(body_fx(-grid % n_bnd_cells:grid % n_cells))
            allocate(body_fy(-grid % n_bnd_cells:grid % n_cells))
            allocate(body_fz(-grid % n_bnd_cells:grid % n_cells))
          end if

          body_fx = 0.0; body_fy = 0.0; body_fz = 0.0

          do s = 1, grid % n_faces
            c1 = grid % faces_c(1,s)
            c2 = grid % faces_c(2,s)

            dotprod = ( (grid % xf(s) - grid % xc(c1)) * grav_x    &
                      + (grid % yf(s) - grid % yc(c1)) * grav_y    &
                      + (grid % zf(s) - grid % zc(c1)) * grav_z )

            body_fx(c1) = body_fx(c1) + flow % density_f(s) * si(s) * dotprod
            body_fy(c1) = body_fy(c1) + flow % density_f(s) * sj(s) * dotprod
            body_fz(c1) = body_fz(c1) + flow % density_f(s) * sk(s) * dotprod

            if (c2 > 0) then
              dotprod = ( (grid % xf(s) - grid % xc(c2)) * grav_x    &
                        + (grid % yf(s) - grid % yc(c2)) * grav_y    &
                        + (grid % zf(s) - grid % zc(c2)) * grav_z )

              body_fx(c2) = body_fx(c2) - flow % density_f(s) * si(s) * dotprod
              body_fy(c2) = body_fy(c2) - flow % density_f(s) * sj(s) * dotprod
              body_fz(c2) = body_fz(c2) - flow % density_f(s) * sk(s) * dotprod

            else

            end if
          end do

        call Comm_Mod_Exchange_Real(grid, body_fx)
        call Comm_Mod_Exchange_Real(grid, body_fy)
        call Comm_Mod_Exchange_Real(grid, body_fz)

        do c = 1, grid % n_cells
          b(c) = b(c) + body_fx(c)
        end do

      case(2)

        do c = 1, grid % n_cells
          b(c) = b(c) + body_fy(c)
        end do

      case(3)

        do c = 1, grid % n_cells
          b(c) = b(c) + body_fz(c)
        end do

    end select

  end if

  end subroutine
