!==============================================================================!
  subroutine Grid_Mod_Get_C1_And_C2_At_Level(grid, level, s, c1, c2)
!------------------------------------------------------------------------------!
!  Get cells surrounding a face (c1 and c2) for a grid, level and face         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: level, s, c1, c2
!==============================================================================!

  ! Cells c1 and c2 on the finest mg level
  c1 = grid % faces_c(1, s)
  c2 = grid % faces_c(2, s)

  ! If not boundary take them from the given coarse level
  if( level > 0 ) then
    if( c2 > 0 ) then
      c1 = grid % level(level) % cell(c1)
      c2 = grid % level(level) % cell(c2)
    end if
  end if

  end subroutine
