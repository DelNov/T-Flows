!==============================================================================!
  subroutine Find_Nearest_Wall_Cell(turb)
!------------------------------------------------------------------------------!
!   The subroutine links interior cells to the closes wall cell. This is       !
!   needed for Standard Smagorinsky SGS model used in 'LES'.                   !
!                                                                              !
!   What if the nearest wall cell is in another processor?                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: HUGE
  use Turb_Mod,  only: Turb_Type
  use Comm_Mod,  only: this_proc
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type) :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: k,  c, nearest_cell
  real                     :: new_distance, old_distance
!==============================================================================!

  grid => turb % pnt_grid

  if(this_proc  < 2)  &
    print *, '# Searching for corresponding wall cells!'

  nearest_cell = 0
  turb % nearest_wall_cell = 0
  old_distance = HUGE
  do c = 1, grid % n_cells
    old_distance = HUGE
    do k = 1, grid % n_cells
      if(grid % cell_near_wall(k)) then
        new_distance = Math_Mod_Distance(                                  &
                                grid % xc(k), grid % yc(k), grid % zc(k),  &
                                grid % xc(c), grid % yc(c), grid % zc(c))
        if(new_distance <= old_distance) then
          nearest_cell =  k
          old_distance = new_distance
        end if 
      end if
    end do
    turb % nearest_wall_cell(c) = nearest_cell
  end do

  if(this_proc < 2) print *, '# Searching finished'

  end subroutine

