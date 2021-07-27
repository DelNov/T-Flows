!==============================================================================!
  subroutine Correct_Fluxes_With_Body_Forces(Flow, Sol)
!------------------------------------------------------------------------------!
!   Calculates body forces (Only due to buoyancy for the time being)           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: t_face_delta => r_face_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: t
  type(Matrix_Type), pointer :: M               ! momentum matrix
  real, contiguous,  pointer :: b(:)
  integer                    :: c1, c2, s
  real                       :: xc1, yc1, zc1, xc2, yc2, zc2
  real                       :: gravity_source, dotprod, a12, dens_h
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  t    => Flow % t
  M    => Sol % M
  b    => Sol % b % val

  !-------------------------------!
  !   For Boussinesq hypothesis   !
  !-------------------------------!
  if(boussinesq) then
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      t_face_delta(s) = t % n(c1) * Grid % f(s)          &
                      + t % n(c2) * (1.0 - Grid % f(s))
      t_face_delta(s) = Flow % t_ref - t_face_delta(s)
    end do
  else
    t_face_delta(1:Grid % n_faces) = 1.0
    Flow % beta                    = 1.0  ! also default from control file
  end if

  ! Correct for Gravity
  if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) >= TINY) then

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      if(c2 > 0) then
        xc1 = Grid % xc(c1)
        yc1 = Grid % yc(c1)
        zc1 = Grid % zc(c1)
        xc2 = Grid % xc(c1) + Grid % dx(s)
        yc2 = Grid % yc(c1) + Grid % dy(s)
        zc2 = Grid % zc(c1) + Grid % dz(s)

        dens_h = 2.0 / ( 1.0 / Flow % density(c1) + 1.0 / Flow % density(c2) )

        ! Interpolate gradients (equation 3.61 in Denner's thesis)
        ! Unit for dotprod: [kg/m/s^2]
        dotprod = 0.5 * dens_h * ( Flow % cell_fx(c1) / Grid % vol(c1)    &
                                 / Flow % density(c1)                     &
                                 + Flow % cell_fx(c2) / Grid % vol(c2)    &
                                 / Flow % density(c2) )                   &
                                 * Grid % dx(s)                           &
                + 0.5 * dens_h * ( Flow % cell_fy(c1) / Grid % vol(c1)    &
                                 / Flow % density(c1)                     &
                                 + Flow % cell_fy(c2) / Grid % vol(c2)    &
                                 / Flow % density(c2) )                   &
                                 * Grid % dy(s)                           &
                + 0.5 * dens_h * ( Flow % cell_fz(c1) / Grid % vol(c1)    &
                                 / Flow % density(c1)                     &
                                 + Flow % cell_fz(c2) / Grid % vol(c2)    &
                                 / Flow % density(c2) )                   &
                                 * Grid % dz(s)

        ! Unit for gravity_source: [kg/m/s^2]
        gravity_source =  &
          (  (Grid % xf(s) - xc1) * grav_x * Flow % beta * t_face_delta(s)   &
           + (Grid % yf(s) - yc1) * grav_y * Flow % beta * t_face_delta(s)   &
           + (Grid % zf(s) - zc1) * grav_z * Flow % beta * t_face_delta(s))  &
             * (Flow % density(c1) - Flow % dens_ref)                        &
         -(  (Grid % xf(s) - xc2) * grav_x * Flow % beta * t_face_delta(s)   &
           + (Grid % yf(s) - yc2) * grav_y * Flow % beta * t_face_delta(s)   &
           + (Grid % zf(s) - zc2) * grav_z * Flow % beta * t_face_delta(s))  &
             * (Flow % density(c2) - Flow % dens_ref)

        ! Units for a12: [m^4s/kg]
        a12 = 0.5 * (  Grid % vol(c1) / M % sav(c1)     &
                     + Grid % vol(c2) / M % sav(c2) ) * M % fc(s)

        ! Unit for gravity_source again: [m^3/s]
        gravity_source =  a12 * (gravity_source - dotprod)

        Flow % v_flux % n(s) = Flow % v_flux % n(s) + gravity_source

        b(c1) = b(c1) - gravity_source
        b(c2) = b(c2) + gravity_source

      end if  ! c2 > 0

    end do

  end if

  end subroutine
