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
  real    :: bulk_flux_x, bulk_xp, bulk_flux_y, bulk_yp, bulk_flux_z, bulk_zp
  integer :: c1, c2, s, i_cel
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, wgt
!==============================================================================!

  ! Set local variables which will not confuse OpenACC
  bulk_flux_x = Flow % bulk % flux_x
  bulk_flux_y = Flow % bulk % flux_y
  bulk_flux_z = Flow % bulk % flux_z
  bulk_xp     = Flow % bulk % xp
  bulk_yp     = Flow % bulk % yp
  bulk_zp     = Flow % bulk % zp

  ! Initialize fluxes to zero
  bulk_flux_x = 0.0
  bulk_flux_y = 0.0
  bulk_flux_z = 0.0

  !----------------------------------------------------------------------!
  !   Summ up volume fluxes [m^3/s] over all faces at monitoring plane   !
  !----------------------------------------------------------------------!

  !$acc parallel loop
  do c1 = 1, grid_n_cells - grid_n_buff_cells
    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)

      if(c1 .lt. c2) then  ! avoid duplicate entries

        !-------!
        !   X   !
        !-------!
        if(grid_xc(c1) <= bulk_xp .and. grid_xc(c1) + grid_dx(s) > bulk_xp) then
          !$acc atomic update
          bulk_flux_x = bulk_flux_x + v_flux(s)
        end if
        if(grid_xc(c1) >= bulk_xp .and. grid_xc(c2) + grid_dx(s) < bulk_xp) then
          !$acc atomic update
          bulk_flux_x = bulk_flux_x - v_flux(s)
        end if

        !-------!
        !   Y   !
        !-------!
        if(grid_yc(c1) <= bulk_yp .and. grid_yc(c1) + grid_dy(s) > bulk_yp) then
          !$acc atomic update
          bulk_flux_y = bulk_flux_y + v_flux(s)
        end if
        if(grid_yc(c1) >= bulk_yp .and. grid_yc(c2) + grid_dy(s) < bulk_yp) then
          !$acc atomic update
          bulk_flux_y = bulk_flux_y - v_flux(s)
        end if

        !-------!
        !   Z   !
        !-------!
        if(grid_zc(c1) <= bulk_zp .and. grid_zc(c1) + grid_dz(s) > bulk_zp) then
          !$acc atomic update
          bulk_flux_z = bulk_flux_z + v_flux(s)
        end if
        if(grid_zc(c1) >= bulk_zp .and. grid_zc(c2) + grid_dz(s) < bulk_zp) then
          !$acc atomic update
          bulk_flux_z = bulk_flux_z - v_flux(s)
        end if

      end if
    end do
    !$acc end loop
  end do
  !$acc end parallel

  call Global % Sum_Real(bulk_flux_x)
  call Global % Sum_Real(bulk_flux_y)
  call Global % Sum_Real(bulk_flux_z)

  ! Set local variables which will not confuse OpenACC
  Flow % bulk % flux_x = bulk_flux_x
  Flow % bulk % flux_y = bulk_flux_y
  Flow % bulk % flux_z = bulk_flux_z

  ! Bulk velocities.  Units: [m^3/s] / [m^2] = [m/s]
  Flow % bulk % u = Flow % bulk % flux_x / (Flow % bulk % area_x + TINY)
  Flow % bulk % v = Flow % bulk % flux_y / (Flow % bulk % area_y + TINY)
  Flow % bulk % w = Flow % bulk % flux_z / (Flow % bulk % area_z + TINY)

  end subroutine
