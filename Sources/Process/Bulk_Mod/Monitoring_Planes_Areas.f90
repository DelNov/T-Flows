!==============================================================================!
  subroutine Bulk_Mod_Monitoring_Planes_Areas(bulk, Grid)
!------------------------------------------------------------------------------!
!   Calculate total surface of the monitoring plane                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Bulk_Type) :: bulk
  type(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, ax_t, ay_t, az_t, wgt
!==============================================================================!

  bulk % area_x = 0.0
  bulk % area_y = 0.0
  bulk % area_z = 0.0

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

      ax_t = abs(Grid % sx(s))
      ay_t = abs(Grid % sy(s))
      az_t = abs(Grid % sz(s))

      ! If the flux is across a buffer face, it is summed up twice.  
      ! The variable "wgt" is here to take care of that.
      wgt = 1.0
      if(c2 > Grid % n_cells - Grid % Comm % n_buff_cells) wgt = 0.5

      !-------!
      !   X   !
      !-------!
      if( (xc1 <= bulk % xp).and.(xc2 > bulk % xp) .or.    &
          (xc2 <= bulk % xp).and.(xc1 > bulk % xp) ) then

        bulk % area_x = bulk % area_x + wgt * ax_t
      end if

      !-------!
      !   Y   !
      !-------!
      if( (yc1 <= bulk % yp).and.(yc2 > bulk % yp) .or. & 
          (yc2 <= bulk % yp).and.(yc1 > bulk % yp) ) then

        bulk % area_y = bulk % area_y + wgt * ay_t
      end if

      !-------!
      !   Z   !
      !-------!
      if( (zc1 <= bulk % zp).and.(zc2 > bulk % zp) .or. & 
          (zc2 <= bulk % zp).and.(zc1 > bulk % zp) ) then

        bulk % area_z = bulk % area_z + wgt * az_t
      end if

    end if

  end do

  call Comm_Mod_Global_Sum_Real(bulk % area_x)
  call Comm_Mod_Global_Sum_Real(bulk % area_y)
  call Comm_Mod_Global_Sum_Real(bulk % area_z)

  end subroutine
