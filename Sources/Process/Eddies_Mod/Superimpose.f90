!==============================================================================!
  subroutine Eddies_Mod_Superimpose(eddies)
!------------------------------------------------------------------------------!
!   Superimpose eddies on velocity field                                     !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  real                      :: xe, ye, ze, xc, yc, zc, uc, vc, wc, re, le, sgn
  real                      :: sig_x, sig_y, sig_z, delta_x, delta_y, delta_z
  integer                   :: e, c, pass
!==============================================================================!

  ! Take aliases to object particle flow around
  grid => eddies % pnt_grid
  flow => eddies % pnt_flow

  do c = -grid % n_bnd_cells, -1
    if(Grid_Mod_Bnd_Cond_Name(grid, c) .eq. eddies % bc_name) then
      flow % u % n(c) = flow % u % b(c)
      flow % v % n(c) = flow % v % b(c)
      flow % w % n(c) = flow % w % b(c)
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

    if( (pass .eq. 1 .or. pass .eq. 2) .and. grid % per_x < NANO) cycle
    if( (pass .eq. 3 .or. pass .eq. 4) .and. grid % per_y < NANO) cycle
    if( (pass .eq. 5 .or. pass .eq. 6) .and. grid % per_z < NANO) cycle

    if( pass .eq. 1 ) delta_x = -grid % per_x
    if( pass .eq. 2 ) delta_x = +grid % per_x
    if( pass .eq. 3 ) delta_y = -grid % per_y
    if( pass .eq. 4 ) delta_y = +grid % per_y
    if( pass .eq. 5 ) delta_z = -grid % per_z
    if( pass .eq. 6 ) delta_z = +grid % per_z

    ! Select cell for each cell randomly
    do c = -grid % n_bnd_cells, -1
      if(Grid_Mod_Bnd_Cond_Name(grid, c) .eq. eddies % bc_name) then

        xc = grid % xc(c)
        yc = grid % yc(c)
        zc = grid % zc(c)

        do e = 1, eddies % n_eddies
          xe  = eddies % eddy(e) % x + delta_x
          ye  = eddies % eddy(e) % y + delta_y
          ze  = eddies % eddy(e) % z + delta_z
          re  = eddies % eddy(e) % radius
          le  = eddies % eddy(e) % length
          sgn = eddies % eddy(e) % sense

          ! Sigmas for Gaussian distribution
          sig_x = re / 2.0
          sig_y = re / 2.0
          sig_z = re / 2.0
          if(eddies % x_plane < HUGE) sig_x = le / 2.0
          if(eddies % y_plane < HUGE) sig_y = le / 2.0
          if(eddies % z_plane < HUGE) sig_z = le / 2.0

          ! Linear rotation
          if(eddies % x_plane < HUGE) then
            uc = 0.0
            vc = sgn * (ze - zc) / re
            wc = sgn * (yc - ye) / re
          else if(eddies % y_plane < HUGE) then
            uc = sgn * (ze - zc) / re
            vc = 0.0
            wc = sgn * (xc - xe) / re
          else if(eddies % z_plane < HUGE) then
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
          uc = uc * ONE_THIRD * 10.0
          vc = vc * ONE_THIRD * 10.0
          wc = wc * ONE_THIRD * 10.0

          flow % u % n(c) = flow % u % n(c) + uc
          flow % v % n(c) = flow % v % n(c) + vc
          flow % w % n(c) = flow % w % n(c) + wc

        end do

      end if
    end do

  end do

  end subroutine
