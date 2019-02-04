!==============================================================================!
  subroutine Solver_Mod_Create(sol, grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type)        :: sol
  type(Grid_Type),  target :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: lev
!==============================================================================!

  sol % pnt_grid => grid

  if(this_proc < 2) print *, '# Determining matrix topology.'

  call Matrix_Mod_Create  (sol % a, grid)
  call Matrix_Mod_Create  (sol % d, grid)
  call Vector_Mod_Allocate(sol % b, grid)

  allocate(sol % a_lev(grid % n_levels))
  allocate(sol % d_lev(grid % n_levels))
  allocate(sol % b_lev(grid % n_levels))
  allocate(sol % x_lev(grid % n_levels))
  allocate(sol % p_lev(grid % n_levels))
  allocate(sol % r_lev(grid % n_levels))

  if(this_proc < 2 .and. grid % n_levels > 1) then
    print '(a10,i3,a22)', ' # Now for', grid % n_levels, ' coarser grid levels.'
  end if

  do lev = 1, grid % n_levels
    call Matrix_Mod_Create_Level  (sol % a_lev(lev), grid, lev)
    call Matrix_Mod_Create_Level  (sol % d_lev(lev), grid, lev)
    call Vector_Mod_Allocate_Level(sol % b_lev(lev), grid, lev)
    call Vector_Mod_Allocate_Level(sol % x_lev(lev), grid, lev)
    call Vector_Mod_Allocate_Level(sol % p_lev(lev), grid, lev)
    call Vector_Mod_Allocate_Level(sol % r_lev(lev), grid, lev)
  end do

  if(this_proc < 2) print *, '# Finished !'

  end subroutine
