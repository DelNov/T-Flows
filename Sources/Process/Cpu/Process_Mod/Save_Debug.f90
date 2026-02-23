!==============================================================================!
  subroutine Save_Debug(Process, Flow, Sol, switch, app)
!------------------------------------------------------------------------------!
!>  Frequently used debug savings.  Maybe it is better to put them all in
!>  one central place, rather than having them scattered around the code.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process
  type(Field_Type),    target :: Flow     !! flow object
  type(Solver_Type)           :: Sol      !! solver object
  character(*)                :: switch   !! what to save?
  character(*)                :: app      !! appending to file name
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  real,    contiguous, pointer :: buff(:)
  integer                      :: c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  Grid => Flow % pnt_grid

  call Work % Connect_Real_Cell(buff)

  if(switch .eq. 'A') then
    do c = Cells_In_Domain()
      buff(c) = Sol % Nat % A % val(Sol % Nat % A % dia(c))
    end do
    call Grid % Save_Debug_Vtu(append = app,        &
                               scalar_cell = buff,  &
                               scalar_name = "A")

  else if(switch .eq. 'B') then
    call Grid % Save_Debug_Vtu(append = app,                       &
                               inside_cell = Sol % Nat % b % val,  &
                               inside_name = "b")

  else if(switch .eq. 'P') then
    call Grid % Save_Debug_Vtu(append = app//"-inside",     &
                               scalar_cell = Flow % p % n,  &
                               scalar_name = "p",           &
                               plot_inside = .true.)
    call Grid % Save_Debug_Vtu(append = app//"-bnd",        &
                               scalar_cell = Flow % p % n,  &
                               scalar_name = "pp",          &
                               plot_inside = .false.)

  else if(switch .eq. 'PP') then
    call Grid % Save_Debug_Vtu(append = app//"-inside",      &
                               scalar_cell = Flow % pp % n,  &
                               scalar_name = "pp",           &
                               plot_inside = .true.)
    call Grid % Save_Debug_Vtu(append = app//"-bnd",         &
                               scalar_cell = Flow % pp % n,  &
                               scalar_name = "pp",           &
                               plot_inside = .false.)

  else if(switch .eq. 'UVW') then
    call Grid % Save_Debug_Vtu(append = app//"-inside",         &
                               vector_cell = (/Flow % u % n,    &
                                               Flow % v % n,    &
                                               Flow % w % n/),  &
                               vector_name = "velocity",        &
                               plot_inside = .true.)
    call Grid % Save_Debug_Vtu(append = app//"-bnd",            &
                               vector_cell = (/Flow % u % n,    &
                                               Flow % v % n,    &
                                               Flow % w % n/),  &
                               vector_name = "velocity",        &
                               plot_inside = .false.)
  end if

  call Work % Disconnect_Real_Cell(buff)

  end subroutine

