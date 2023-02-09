!==============================================================================!
  subroutine Calculate_Fluxes(Flow, v_flux)
!------------------------------------------------------------------------------!
!   Calculate mass fluxes through whole domain.                                  !
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
  real                     :: dens_f
  real                     :: dens_are_x, dens_are_y, dens_are_z  ! [kg/m]
  real                     :: dens_avg_x, dens_avg_y, dens_avg_z  ! [kg/m^3]
!==============================================================================!

  ! Take some aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk

  bulk % flux_x = 0.0
  bulk % flux_y = 0.0
  bulk % flux_z = 0.0

  dens_are_x = 0.0
  dens_are_y = 0.0
  dens_are_z = 0.0

  !-------------------------------------------------------------------------!
  !   Summ up mass fluxes [kg/s] over all faces at monitoring plane         !
  !   (Resulting mass flux will be for the whole domain, still in [kg/s])   !
  !   In addition, sum up face densities time faca areas [kg/m]             !
  !-------------------------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Density at the face
    dens_f = Flow % density(c1) *        Grid % fw(s)  &
           + Flow % density(c2) * (1.0 - Grid % fw(s))

    if(Grid % Comm % cell_proc(c1) .eq. this_proc) then

      xc1 = Grid % xc(c1)
      yc1 = Grid % yc(c1)
      zc1 = Grid % zc(c1)
      xc2 = Grid % xc(c1) + Grid % dx(s)
      yc2 = Grid % yc(c1) + Grid % dy(s)
      zc2 = Grid % zc(c1) + Grid % dz(s)

      ! If the flux is across a buffer face, it is summed up twice.  
      ! The variable "wgt" is here to take care of that.
      wgt = 1.0
      if(c2 > Grid % n_cells - Grid % Comm % n_buff_cells) wgt = 0.5

      !-------!
      !   X   !
      !-------!
      if((xc1 <= bulk % xp).and.(xc2 > bulk % xp)) then
        bulk % flux_x = bulk % flux_x + wgt * v_flux(s) * dens_f
        dens_are_x = dens_are_x + wgt * dens_f * abs(Grid % sx(s))
      end if
      if((xc2 < bulk % xp).and.(xc1 >= bulk % xp)) then
        bulk % flux_x = bulk % flux_x - wgt * v_flux(s) * dens_f
        dens_are_x = dens_are_x + wgt * dens_f * abs(Grid % sx(s))
      end if

      !-------!
      !   Y   !
      !-------!
      if((yc1 <= bulk % yp).and.(yc2 > bulk % yp)) then
        bulk % flux_y = bulk % flux_y + wgt * v_flux(s) * dens_f
        dens_are_y = dens_are_y + wgt * dens_f * abs(Grid % sy(s))
      end if
      if((yc2 < bulk % yp).and.(yc1 >= bulk % yp)) then
        bulk % flux_y = bulk % flux_y - wgt * v_flux(s) * dens_f
        dens_are_y = dens_are_y + wgt * dens_f * abs(Grid % sy(s))
      end if

      !-------!
      !   Z   !
      !-------!
      if((zc1 <= bulk % zp).and.(zc2 > bulk % zp)) then
        bulk % flux_z = bulk % flux_z + wgt * v_flux(s) * dens_f
        dens_are_z = dens_are_z + wgt * dens_f * abs(Grid % sz(s))
      end if
      if((zc2 < bulk % zp).and.(zc1 >= bulk % zp)) then
        bulk % flux_z = bulk % flux_z - wgt * v_flux(s) * dens_f
        dens_are_z = dens_are_z + wgt * dens_f * abs(Grid % sz(s))
      end if

    end if
  end do

  ! Calculate average densities in three characterstic planes
  call Comm_Mod_Global_Sum_Real(dens_are_x)
  call Comm_Mod_Global_Sum_Real(dens_are_y)
  call Comm_Mod_Global_Sum_Real(dens_are_z)
  dens_avg_x = dens_are_x / (bulk % area_x + NANO)
  dens_avg_y = dens_are_y / (bulk % area_y + NANO)
  dens_avg_z = dens_are_z / (bulk % area_z + NANO)

  call Comm_Mod_Global_Sum_Real(bulk % flux_x)
  call Comm_Mod_Global_Sum_Real(bulk % flux_y)
  call Comm_Mod_Global_Sum_Real(bulk % flux_z)

  ! Bulk velocities.  Units: [kg/s] / [kg/m^3] / [m^2] = [m/s]
  bulk % u = bulk % flux_x / (dens_avg_x * bulk % area_x + TINY)
  bulk % v = bulk % flux_y / (dens_avg_y * bulk % area_y + TINY)
  bulk % w = bulk % flux_z / (dens_avg_z * bulk % area_z + TINY)

  end subroutine
