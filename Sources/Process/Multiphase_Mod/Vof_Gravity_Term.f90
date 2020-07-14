!==============================================================================!
  subroutine Multiphase_Mod_Vof_Gravity_Term(mult)
!------------------------------------------------------------------------------!
!   Computes Gravity sources for Momentum Equation                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  integer                   :: s, c, c1, c2
  real                      :: dotprod, epsloc
  real,             pointer :: si(:), sj(:), sk(:)
  real                      :: corr_x, corr_y, corr_z
!==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => mult % pnt_grid
  vof    => mult % vof
  si     => grid % sx
  sj     => grid % sy
  sk     => grid % sz

  epsloc = epsilon(epsloc)

  if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

    if( .not. allocated(mult % body_fx)) then
      allocate(mult % body_fx(-grid % n_bnd_cells:grid % n_cells))
      allocate(mult % body_fy(-grid % n_bnd_cells:grid % n_cells))
      allocate(mult % body_fz(-grid % n_bnd_cells:grid % n_cells))
    end if

    mult % body_fx = 0.0; mult % body_fy = 0.0; mult % body_fz = 0.0

    ! At boundaries
    do s = 1, grid % n_bnd_faces

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
    end do

    ! At interior faces
    do s = grid % n_bnd_faces + 1, grid % n_faces

      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      dotprod = ( (grid % xf(s) - grid % xc(c1)) * grav_x    &
                + (grid % yf(s) - grid % yc(c1)) * grav_y    &
                + (grid % zf(s) - grid % zc(c1)) * grav_z )

      mult % body_fx(c1) = mult % body_fx(c1)                      &
                         + flow % density_f(s) * si(s) * dotprod
      mult % body_fy(c1) = mult % body_fy(c1)                      &
                         + flow % density_f(s) * sj(s) * dotprod
      mult % body_fz(c1) = mult % body_fz(c1)                      &
                         + flow % density_f(s) * sk(s) * dotprod

      ! Correction for periodic faces:
      call Grid_Mod_Correction_Periodicity(grid, s,   &
                                           corr_x, corr_y, corr_z)

      dotprod = ( ( grid % xf(s) + corr_x - grid % xc(c2)) * grav_x     &
                + ( grid % yf(s) + corr_y - grid % yc(c2)) * grav_y     &
                + ( grid % zf(s) + corr_z - grid % zc(c2)) * grav_z )

      mult % body_fx(c2) = mult % body_fx(c2)                      &
                         - flow % density_f(s) * si(s) * dotprod
      mult % body_fy(c2) = mult % body_fy(c2)                      &
                         - flow % density_f(s) * sj(s) * dotprod
      mult % body_fz(c2) = mult % body_fz(c2)                      &
                         - flow % density_f(s) * sk(s) * dotprod

    end do

    call Grid_Mod_Exchange_Cells_Real(grid, mult % body_fx)
    call Grid_Mod_Exchange_Cells_Real(grid, mult % body_fy)
    call Grid_Mod_Exchange_Cells_Real(grid, mult % body_fz)

  end if

  end subroutine
