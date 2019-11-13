!==============================================================================!
  subroutine Eddies_Mod_Superimpose(eddies)
!------------------------------------------------------------------------------!
!   Create memory for eddies                                                 !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  real                      :: xe, ye, ze, xc, yc, zc, uc, vc, wc, re, le, sgn
  real                      :: sig_x, sig_y, sig_z
  integer                   :: e, c
!==============================================================================!

  ! Take aliases to object particle flow around
  grid => eddies % pnt_grid
  flow => eddies % pnt_flow

  ! Select cell for each cell randomly
  do c = -grid % n_bnd_cells, -1
    if(Grid_Mod_Bnd_Cond_Name(grid, c) .eq. eddies % bc_name) then

      xc = grid % xc(c)
      yc = grid % yc(c)
      zc = grid % zc(c)

      do e = 1, eddies % n_eddies
        xe  = eddies % eddy(e) % x
        ye  = eddies % eddy(e) % y
        ze  = eddies % eddy(e) % z
        re  = eddies % eddy(e) % radius
        le  = eddies % eddy(e) % length
        sgn = eddies % eddy(e) % sgn

        ! Sigmas for Gaussian distribution
        sig_x = le / 2.0
        sig_y = re / 2.0
        sig_z = re / 2.0

        ! Velocity components assuming flow in x direction

        ! Linear rotation
        uc = 0.0
        vc = sgn * (ze - zc) / re
        wc = sgn * (yc - ye) / re

        ! Damp rotation with Gaussian (like) distribution
        vc = vc * exp(-0.5*((xc-xe)/sig_x)**2)
        vc = vc * exp(-0.5*((yc-ye)/sig_y)**2)
        vc = vc * exp(-0.5*((zc-ze)/sig_z)**2)

        wc = wc * exp(-0.5*((xc-xe)/sig_x)**2)
        wc = wc * exp(-0.5*((yc-ye)/sig_y)**2)
        wc = wc * exp(-0.5*((zc-ze)/sig_z)**2)

        ! Scale them so that intensity turns to be one
        vc = vc * ONE_THIRD * 10.0
        wc = wc * ONE_THIRD * 10.0

        flow % v % n(c) = flow % v % n(c) + vc
        flow % w % n(c) = flow % w % n(c) + wc

      end do

    end if
  end do

  end subroutine
