!==============================================================================!
  subroutine Superimpose_Eddies(Eddies)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to superimpose the effects of eddies onto the
!>  velocity field, useful for some aspects scale-resolving simulations.
!>  The subroutine effectively superimposes rotational velocities induced by
!>  the eddies onto the existing velocity field at the boundary cells.
!>  This superimposition simulates the effect of turbulence entering the
!>  computational domain.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine first resets the velocity field at the boundary cells     !
!     to their baseline values.                                                !
!   * It then goes through different passes, each corresponding to a p         !
!     otential periodic direction in the computational domain. In each pass,   !
!     the subroutine adjusts the positions of the eddies based on the          !
!     periodicity of the grid. This is done to manage eddies passing through   !
!     periodic boundary conditions in the domain,                              !
!   * For each boundary cell within the specified boundary condition plane,    !
!     the subroutine calculates the contributions from all eddies. It          !
!     considers each eddy's position, radius, length, and sense of rotation.   !
!   * The influence of each eddy on a boundary cell is calculated using a      !
!     Gaussian-like distribution to dampen the velocity induced by the eddy    !
!     based on its distance from the cell. The Gaussian distribution ensures   !
!     that the influence of each eddy diminishes with distance.                !
!   * The velocities induced by the eddies are then scaled based on the        !
!     intensity parameter of the eddies structure and added to the existing    !
!     velocities of the boundary cells.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Eddies_Type), target :: Eddies  !! parent class; collection of eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  real                      :: xe, ye, ze, xc, yc, zc, uc, vc, wc, re, le, sgn
  real                      :: sig_x, sig_y, sig_z, delta_x, delta_y, delta_z
  integer                   :: e, c, pass
!==============================================================================!

  ! Take aliases to object particle Flow around
  Grid => Eddies % pnt_grid
  Flow => Eddies % pnt_flow

  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name_At_Cell(c) .eq. Eddies % bc_name) then
      Flow % u % n(c) = Flow % u % b(c)
      Flow % v % n(c) = Flow % v % b(c)
      Flow % w % n(c) = Flow % w % b(c)
    end if
  end do

  !----------------------------------------------------------!
  !     Browse through all possible periodic directions      !
  !   (for each direction, it will shift eddies positions)   !
  !----------------------------------------------------------!
  do pass = 1, 7
    delta_x = 0
    delta_y = 0
    delta_z = 0

    if( (pass .eq. 1 .or. pass .eq. 2) .and. Grid % per_x < NANO) cycle
    if( (pass .eq. 3 .or. pass .eq. 4) .and. Grid % per_y < NANO) cycle
    if( (pass .eq. 5 .or. pass .eq. 6) .and. Grid % per_z < NANO) cycle

    if( pass .eq. 1 ) delta_x = -Grid % per_x
    if( pass .eq. 2 ) delta_x = +Grid % per_x
    if( pass .eq. 3 ) delta_y = -Grid % per_y
    if( pass .eq. 4 ) delta_y = +Grid % per_y
    if( pass .eq. 5 ) delta_z = -Grid % per_z
    if( pass .eq. 6 ) delta_z = +Grid % per_z

    ! Select cell for each cell randomly
    do c = -Grid % n_bnd_cells, -1
      if(Grid % Bnd_Cond_Name_At_Cell(c) .eq. Eddies % bc_name) then

        xc = Grid % xc(c)
        yc = Grid % yc(c)
        zc = Grid % zc(c)

        do e = 1, Eddies % n_eddies
          xe  = Eddies % eddy(e) % x + delta_x
          ye  = Eddies % eddy(e) % y + delta_y
          ze  = Eddies % eddy(e) % z + delta_z
          re  = Eddies % eddy(e) % radius
          le  = Eddies % eddy(e) % length
          sgn = Eddies % eddy(e) % sense

          ! Sigmas for Gaussian distribution
          sig_x = re / 2.0
          sig_y = re / 2.0
          sig_z = re / 2.0
          if(Eddies % x_plane < HUGE) sig_x = le / 2.0
          if(Eddies % y_plane < HUGE) sig_y = le / 2.0
          if(Eddies % z_plane < HUGE) sig_z = le / 2.0

          ! Linear rotation
          if(Eddies % x_plane < HUGE) then
            uc = 0.0
            vc = sgn * (ze - zc) / re
            wc = sgn * (yc - ye) / re
          else if(Eddies % y_plane < HUGE) then
            uc = sgn * (ze - zc) / re
            vc = 0.0
            wc = sgn * (xc - xe) / re
          else if(Eddies % z_plane < HUGE) then
            uc = sgn * (yc - ye) / re
            vc = sgn * (xe - xc) / re
            wc = 0.0
          end if

          ! Damp rotation with Gaussian (like) distribution
          uc = uc * exp(-0.5*((xc-xe)/sig_x)**2)
          uc = uc * exp(-0.5*((yc-ye)/sig_y)**2)
          uc = uc * exp(-0.5*((zc-ze)/sig_z)**2)

          vc = vc * exp(-0.5*((xc-xe)/sig_x)**2)
          vc = vc * exp(-0.5*((yc-ye)/sig_y)**2)
          vc = vc * exp(-0.5*((zc-ze)/sig_z)**2)

          wc = wc * exp(-0.5*((xc-xe)/sig_x)**2)
          wc = wc * exp(-0.5*((yc-ye)/sig_y)**2)
          wc = wc * exp(-0.5*((zc-ze)/sig_z)**2)

          ! Scale them so that intensity turns to be equal to one
          uc = uc * ONE_THIRD * 10.0 * Eddies % intensity
          vc = vc * ONE_THIRD * 10.0 * Eddies % intensity
          wc = wc * ONE_THIRD * 10.0 * Eddies % intensity

          Flow % u % n(c) = Flow % u % n(c) + uc
          Flow % v % n(c) = Flow % v % n(c) + vc
          Flow % w % n(c) = Flow % w % n(c) + wc

        end do

      end if
    end do

  end do

  end subroutine
