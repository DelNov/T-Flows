!==============================================================================!
  subroutine Bulk_Mod_Monitoring_Planes_Areas(bulk, grid)
!------------------------------------------------------------------------------!
!   Calculate total surface of the monitoring plane                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Bulk_Type) :: bulk
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, m
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, ax_t, ay_t, az_t, wgt
!==============================================================================!

  bulk % area_x = 0.0
  bulk % area_y = 0.0
  bulk % area_z = 0.0

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

      ax_t = abs(grid % sx(s))
      ay_t = abs(grid % sy(s))
      az_t = abs(grid % sz(s))

      ! If the flux is across a buffer face, it is summed up twice.  
      ! The variable "wgt" is here to take care of that.
      wgt = 1.0
      if(c2 > grid % n_cells - grid % comm % n_buff_cells) wgt = 0.5

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
