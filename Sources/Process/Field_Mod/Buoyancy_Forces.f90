!==============================================================================!
  subroutine Field_Mod_Buoyancy_Forces(flow)
!------------------------------------------------------------------------------!
!   Calculates buoyancy with or without Boussinesq hypothesis.  That means     !
!   it is used for both single phase flows with density assumed constant and   !
!   for multiphase flows with VOF and large density variations.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: dens_f => r_face_01
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                                    (allocates Work_Mod)           !
!     |                                                                        !
!     +----> Compute_Momentum                   (doesn't use Work_Mod)         !
!             |                                                                !
!             +----> Field_Mod_Buoyancy_Forces  (safe to use r_face_01)        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: t
  integer                    :: c, c1, c2, s
  real                       :: xc1, yc1, zc1, xc2, yc2, zc2, dotprod, temp_f
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t

  ! Initialize
  flow % face_fx(:) = 0.0
  flow % face_fy(:) = 0.0
  flow % face_fz(:) = 0.0
  flow % cell_fx(:) = 0.0
  flow % cell_fy(:) = 0.0
  flow % cell_fz(:) = 0.0

  !-------------------------------!
  !   For Boussinesq hypothesis   !
  !-------------------------------!
  if(boussinesq) then
    do s = 1, grid % n_faces
      c1  = grid % faces_c(1, s)
      c2  = grid % faces_c(2, s)

      ! Find temperature and density with linear interpolation.
      ! (One could consider a more accurate interpolation with
      !  gradients, or even just running it through Grad_Gauss
      !  to get iterativelly improved face temperature values)
      temp_f =        grid % f(s)  * t % n(c1)  &
             + (1.0 - grid % f(s)) * t % n(c2)

      ! Density is constant for Boussinesq hypothesis
      ! this linear interpolation should do just fine
      dens_f(s) =        grid % fw(s)  * flow % density(c1)  &
                + (1.0 - grid % fw(s)) * flow % density(c2)

      ! Transform density into one assumed by Boussinesq
      ! Units: kg/m^3 * K * 1/K = kg/m^3
      dens_f(s) = dens_f(s) * (flow % t_ref - temp_f) * flow % beta
    end do

  !---------------------------------------------------------------------!
  !   no Boussinesq hypothesis, this will include multiphase flow and   !
  !    probably very abbrupt changes in density, go for harmonic mean   !
  !---------------------------------------------------------------------!
  else
    do s = 1, grid % n_faces
      c1  = grid % faces_c(1, s)
      c2  = grid % faces_c(2, s)

      if(c2 > 0) then
        dens_f(s) = Math_Mod_Harmonic_Mean(flow % density(c1),  &
                                           flow % density(c2))
      else
        dens_f(s) = flow % density(c1)
      end if

      dens_f(s) = dens_f(s) - flow % dens_ref
    end do
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

      ! Units in next three expressions which follow
      ! kg/m^3 * m/s^2 * m^2 * m = kg m/s^2 = N
      flow % cell_fx(c1) = flow % cell_fx(c1)          &
                         + dens_f(s) * grav_x          &
                         * grid % sx(s) * (grid % xf(s) - xc1)

      flow % cell_fy(c1) = flow % cell_fy(c1)          &
                         + dens_f(s) * grav_y          &
                         * grid % sy(s) * (grid % yf(s) - yc1)

      flow % cell_fz(c1) = flow % cell_fz(c1)          &
                         + dens_f(s) * grav_z          &
                         * grid % sz(s) * (grid % zf(s) - zc1)

      if(c2 > 0) then
        xc2 = grid % xc(c1) + grid % dx(s)
        yc2 = grid % yc(c1) + grid % dy(s)
        zc2 = grid % zc(c1) + grid % dz(s)

        ! Units in next three expressions which follow
        ! kg/m^3 * m/s^2 * m^2 * m = kg m/s^2 = N
        flow % cell_fx(c2) = flow % cell_fx(c2)          &
                           - dens_f(s) * grav_x          &
                           * grid % sx(s) * (grid % xf(s) - xc2)

        flow % cell_fy(c2) = flow % cell_fy(c2)          &
                           - dens_f(s) * grav_y          &
                           * grid % sy(s) * (grid % yf(s) - yc2)

        flow % cell_fz(c2) = flow % cell_fz(c2)          &
                           - dens_f(s) * grav_z          &
                           * grid % sz(s) * (grid % zf(s) - zc2)
      end if
    end do

  end if  ! boussinesq

  call Grid_Mod_Exchange_Cells_Real(grid, flow % cell_fx)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % cell_fy)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % cell_fz)

  end subroutine
