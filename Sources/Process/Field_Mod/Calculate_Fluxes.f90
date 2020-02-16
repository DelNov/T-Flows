!==============================================================================!
  subroutine Field_Mod_Calculate_Fluxes(flow, flux)
!------------------------------------------------------------------------------!
!   Calculate mass fluxes through whole domain.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  real, dimension(:)       :: flux
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  integer                  :: c1, c2, s
  real                     :: xc1, yc1, zc1, xc2, yc2, zc2, wgt
  real                     :: dens_are_x, dens_are_y, dens_are_z  ! [kg/m]
  real                     :: dens_avg_x, dens_avg_y, dens_avg_z  ! [kg/m^3]
!==============================================================================!

  ! Take some aliases
  grid => flow % pnt_grid
  bulk => flow % bulk

  bulk % flux_x = 0.0
  bulk % flux_y = 0.0
  bulk % flux_z = 0.0

  !-------------------------------------------------------------------------!
  !   Summ up mass fluxes [kg/s] over all faces at monitoring plane         !
  !   (Resulting mass flux will be for the whole domain, still in [kg/s])   !
  !   In addition, sum up face densities time faca areas [kg/m]             !
  !-------------------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 > 0) then

      xc1 = grid % xc(c1)
      yc1 = grid % yc(c1)
      zc1 = grid % zc(c1)
      xc2 = grid % xc(c1) + grid % dx(s)
      yc2 = grid % yc(c1) + grid % dy(s)
      zc2 = grid % zc(c1) + grid % dz(s)

      ! If the flux is across a buffer face, it is summed up twice.  
      ! The variable "wgt" is here to take care of that.
      wgt = 1.0
      if(c2 > grid % n_cells - grid % comm % n_buff_cells) wgt = 0.5

      !-------!
      !   X   !
      !-------!
      if((xc1 <= bulk % xp).and.(xc2 > bulk % xp)) then
        bulk % flux_x = bulk % flux_x + wgt * flux(s)
        dens_are_x = dens_are_x + wgt * flow % density_f(s) * abs(grid % sx(s))
      end if
      if((xc2 < bulk % xp).and.(xc1 >= bulk % xp)) then
        bulk % flux_x = bulk % flux_x - wgt * flux(s)
        dens_are_x = dens_are_x + wgt * flow % density_f(s) * abs(grid % sx(s))
      end if

      !-------!
      !   Y   !
      !-------!
      if((yc1 <= bulk % yp).and.(yc2 > bulk % yp)) then
        bulk % flux_y = bulk % flux_y + wgt * flux(s)
        dens_are_y = dens_are_y + wgt * flow % density_f(s) * abs(grid % sy(s))
      end if
      if((yc2 < bulk % yp).and.(yc1 >= bulk % yp)) then
        bulk % flux_y = bulk % flux_y - wgt * flux(s)
        dens_are_y = dens_are_y + wgt * flow % density_f(s) * abs(grid % sy(s))
      end if

      !-------!
      !   Z   !
      !-------!
      if((zc1 <= bulk % zp).and.(zc2 > bulk % zp)) then
        bulk % flux_z = bulk % flux_z + wgt * flux(s)
        dens_are_z = dens_are_z + wgt * flow % density_f(s) * abs(grid % sz(s))
      end if
      if((zc2 < bulk % zp).and.(zc1 >= bulk % zp)) then
        bulk % flux_z = bulk % flux_z - wgt * flux(s)
        dens_are_z = dens_are_z + wgt * flow % density_f(s) * abs(grid % sz(s))
      end if

    end if
  end do

  ! Calculate average densities in three characterstic planes
  call Comm_Mod_Global_Sum_Real(dens_are_x)
  call Comm_Mod_Global_Sum_Real(dens_are_y)
  call Comm_Mod_Global_Sum_Real(dens_are_z)
  dens_avg_x = dens_are_x / bulk % area_x
  dens_avg_y = dens_are_y / bulk % area_y
  dens_avg_z = dens_are_z / bulk % area_z

  call Comm_Mod_Global_Sum_Real(bulk % flux_x)
  call Comm_Mod_Global_Sum_Real(bulk % flux_y)
  call Comm_Mod_Global_Sum_Real(bulk % flux_z)

  ! Bulk velocities.  Units: [kg/s] / [kg/m] = [m/s]
  bulk % u = bulk % flux_x / (dens_avg_x + TINY)
  bulk % v = bulk % flux_y / (dens_avg_y + TINY)
  bulk % w = bulk % flux_z / (dens_avg_z + TINY)

  end subroutine
