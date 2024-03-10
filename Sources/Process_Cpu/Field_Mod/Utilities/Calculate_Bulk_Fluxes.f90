!==============================================================================!
  subroutine Calculate_Bulk_Fluxes(Flow, v_flux)
!------------------------------------------------------------------------------!
!   Calculate volume fluxes through whole domain.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  real, dimension(:)        :: v_flux  ! volume flux
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: c1, c2, s
  real                     :: xc1, yc1, zc1, xc2, yc2, zc2, wgt
!==============================================================================!

  ! Take some aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk

  bulk % flux_x = 0.0
  bulk % flux_y = 0.0
  bulk % flux_z = 0.0

  !----------------------------------------------------------------------!
  !   Summ up volume fluxes [m^3/s] over all faces at monitoring plane   !
  !----------------------------------------------------------------------!
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    xc1 = Grid % xc(c1)
    yc1 = Grid % yc(c1)
    zc1 = Grid % zc(c1)
    xc2 = Grid % xc(c1) + Grid % dx(s)
    yc2 = Grid % yc(c1) + Grid % dy(s)
    zc2 = Grid % zc(c1) + Grid % dz(s)

    ! If the flux is across a buffer face, it is summed up twice.
    ! The variable "wgt" is here to take care of that.
    wgt = 1.0
    if(.not. Cell_In_This_Proc(c2)) wgt = 0.5

    !-------!
    !   X   !
    !-------!
    if((xc1 <= bulk % xp).and.(xc2 > bulk % xp)) then
      bulk % flux_x = bulk % flux_x + wgt * v_flux(s)
    end if
    if((xc2 < bulk % xp).and.(xc1 >= bulk % xp)) then
      bulk % flux_x = bulk % flux_x - wgt * v_flux(s)
    end if

    !-------!
    !   Y   !
    !-------!
    if((yc1 <= bulk % yp).and.(yc2 > bulk % yp)) then
      bulk % flux_y = bulk % flux_y + wgt * v_flux(s)
    end if
    if((yc2 < bulk % yp).and.(yc1 >= bulk % yp)) then
      bulk % flux_y = bulk % flux_y - wgt * v_flux(s)
    end if

    !-------!
    !   Z   !
    !-------!
    if((zc1 <= bulk % zp).and.(zc2 > bulk % zp)) then
      bulk % flux_z = bulk % flux_z + wgt * v_flux(s)
    end if
    if((zc2 < bulk % zp).and.(zc1 >= bulk % zp)) then
      bulk % flux_z = bulk % flux_z - wgt * v_flux(s)
    end if

  end do

  call Global % Sum_Real(bulk % flux_x)
  call Global % Sum_Real(bulk % flux_y)
  call Global % Sum_Real(bulk % flux_z)

  ! Bulk velocities.  Units: [m^3/s] / [m^2] = [m/s]
  bulk % u = bulk % flux_x / (bulk % area_x + TINY)
  bulk % v = bulk % flux_y / (bulk % area_y + TINY)
  bulk % w = bulk % flux_z / (bulk % area_z + TINY)

  end subroutine
