!==============================================================================!
  subroutine Find_Nearest_Wall_Cell(Turb)
!------------------------------------------------------------------------------!
!   The subroutine links interior cells to the closes wall cell. This is       !
!   needed for Standard Smagorinsky SGS model used in 'LES'.                   !
!                                                                              !
!   What if the nearest wall cell is in another processor?                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: min_dis => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c1, c2, s
!==============================================================================!

  Grid => Turb % pnt_grid

  if(this_proc  < 2)  &
    print *, '# Searching for corresponding wall cells!'

  min_dis(:) = HUGE
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 > 0) then

      ! Cell c1 is near wall cell
      if(Grid % cell_near_wall(c1) .and. .not. Grid % cell_near_wall(c2)) then
        if(min_dis(c2) > Grid % d(s)) then
          min_dis(c2) = Grid % d(s)
          Turb % nearest_wall_cell(c2) = c1
        end if
      end if

      ! Cell c2 is near wall cell
      if(Grid % cell_near_wall(c2) .and. .not. Grid % cell_near_wall(c1)) then
        if(min_dis(c1) > Grid % d(s)) then
          min_dis(c1) = Grid % d(s)
          Turb % nearest_wall_cell(c1) = c2
        end if
      end if

    end if
  end do

  if(this_proc < 2) print *, '# Searching finished'

  end subroutine

