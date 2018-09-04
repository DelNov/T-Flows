!==============================================================================!
  subroutine Bulk_Mod_Monitoring_Planes_Areas(grid, bulk)
!------------------------------------------------------------------------------!
!   Calculate total surface of the monitoring plane                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  type(Bulk_Type) :: bulk
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, m
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, ax_t, ay_t, az_t
!==============================================================================!

  bulk % area_x = 0.0
  bulk % area_y = 0.0
  bulk % area_z = 0.0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  < 0) then
      if (Grid_Mod_Bnd_Cond_Type(grid, c2) .ne. BUFFER) then
        cycle ! skip current surface
      end if
    end if

    xc1 = grid % xc(c1)
    yc1 = grid % yc(c1)
    zc1 = grid % zc(c1)
    xc2 = grid % xc(c1) + grid % dx(s)
    yc2 = grid % yc(c1) + grid % dy(s)
    zc2 = grid % zc(c1) + grid % dz(s)

    ax_t = abs(grid % sx(s))
    ay_t = abs(grid % sy(s))
    az_t = abs(grid % sz(s))

    !-------!
    !   X   !
    !-------!
    if( (xc1 <= bulk % xp).and.(xc2 > bulk % xp) .or.    &
        (xc2 <= bulk % xp).and.(xc1 > bulk % xp) ) then

      ! Watch out: buffer cell faces will be counted twice
      if(c2  < 0) then
        if (Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. BUFFER) then
            bulk % area_x = bulk % area_x + 0.5 * ax_t
        end if
      else
        bulk % area_x = bulk % area_x + ax_t    
      end if
    end if

    !-------!
    !   Y   !
    !-------!
    if( (yc1 <= bulk % yp).and.(yc2 > bulk % yp) .or. & 
        (yc2 <= bulk % yp).and.(yc1 > bulk % yp) ) then

      ! Watch out: buffer cell faces will be counted twice
      if(c2  < 0) then
        if (Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. BUFFER) then
          bulk % area_y = bulk % area_y + 0.5 * ay_t
        end if
      else
        bulk % area_y = bulk % area_y + ay_t    
      end if
    end if

    !-------!
    !   Z   !
    !-------!
    if( (zc1 <= bulk % zp).and.(zc2 > bulk % zp) .or. & 
        (zc2 <= bulk % zp).and.(zc1 > bulk % zp) ) then

      ! Watch out: buffer cell faces will be counted twice
      if(c2  < 0) then
        if (Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. BUFFER) then
          bulk % area_z = bulk % area_z + 0.5 * az_t
        end if
      else
        bulk % area_z = bulk % area_z + az_t    
      end if
    end if

  end do

  call Comm_Mod_Global_Sum_Real(bulk % area_x)
  call Comm_Mod_Global_Sum_Real(bulk % area_y)
  call Comm_Mod_Global_Sum_Real(bulk % area_z)

  end subroutine
