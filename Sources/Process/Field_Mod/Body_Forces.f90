!==============================================================================!
  subroutine Field_Mod_Body_Forces(flow)
!------------------------------------------------------------------------------!
!   Calculates body forces (Only due to buoyancy for the time being)           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: t_face_delta => r_face_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: t
  integer                    :: c1, c2, s, c
  real                       :: xc1, yc1, zc1, xc2, yc2, zc2, dotprod, dens_f
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t

  ! Initialize
  flow % body_fx(:) = 0.0
  flow % body_fy(:) = 0.0
  flow % body_fz(:) = 0.0

  !-------------------------------!
  !   For Boussinesq hypothesis   !
  !-------------------------------!
  if(buoyancy) then
    do s = 1, grid % n_faces
      c1  = grid % faces_c(1, s)
      c2  = grid % faces_c(2, s)
      t_face_delta(s) = t % n(c1) * grid % f(s)          &
                      + t % n(c2) * (1.0 - grid % f(s))
      t_face_delta(s) = t_face_delta(s) - t_ref
    end do

  !------------------------------!
  !   no Boussinesq hypothesis   !
  !------------------------------!
  else
    t_face_delta(1:grid % n_faces) = 1.0
    flow % beta                    = 1.0  ! also default from control file
  end if

  !--------------------!
  !   Buoyancy force   !
  !--------------------!
  if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

    ! Calculate
    do s = 1, grid % n_faces
      c1  = grid % faces_c(1, s)
      c2  = grid % faces_c(2, s)
      xc1 = grid % xc(c1)
      yc1 = grid % yc(c1)
      zc1 = grid % zc(c1)

      dens_f = flow % density(c1) *        grid % fw(s)  &
             + flow % density(c2) * (1.0 - grid % fw(s))

      dotprod = (   (grid % xf(s) - xc1) * grav_x    &
                  + (grid % yf(s) - yc1) * grav_y    &
                  + (grid % zf(s) - zc1) * grav_z )

      flow % body_fx(c1) = flow % body_fx(c1)          &
                         + dens_f * t_face_delta(s)    &
                         * flow % beta * grid % sx(s) * dotprod

      flow % body_fy(c1) = flow % body_fy(c1)          &
                         + dens_f * t_face_delta(s)    &
                         * flow % beta * grid % sy(s) * dotprod

      flow % body_fz(c1) = flow % body_fz(c1)          &
                         + dens_f * t_face_delta(s)    &
                         * flow % beta * grid % sz(s) * dotprod

      if(c2 > 0) then
        xc2 = grid % xc(c1) + grid % dx(s)
        yc2 = grid % yc(c1) + grid % dy(s)
        zc2 = grid % zc(c1) + grid % dz(s)
        dotprod = (   (grid % xf(s) - xc2) * grav_x    &
                    + (grid % yf(s) - yc2) * grav_y    &
                    + (grid % zf(s) - zc2) * grav_z )

        flow % body_fx(c2) = flow % body_fx(c2)          &
                           - dens_f * t_face_delta(s)    &
                           * flow % beta * grid % sx(s) * dotprod

        flow % body_fy(c2) = flow % body_fy(c2)          &
                           - dens_f * t_face_delta(s)    &
                           * flow % beta * grid % sy(s) * dotprod

        flow % body_fz(c2) = flow % body_fz(c2)          &
                           - dens_f * t_face_delta(s)    &
                           * flow % beta * grid % sz(s) * dotprod
      end if
    end do

  end if  ! buoyancy

  call Grid_Mod_Exchange_Cells_Real(grid, flow % body_fx)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % body_fy)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % body_fz)

  ! Was used for debugging
  ! do c = 1, grid % n_cells
  !   flow % pot % n(c) = flow % body_fz(c)
  ! end do

  end subroutine
