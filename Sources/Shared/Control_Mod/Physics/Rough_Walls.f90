!==============================================================================!
  subroutine Control_Mod_Rough_Walls(rough_walls, verbose)
!------------------------------------------------------------------------------!
!   Reading wall roughness from the control file.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: rough_walls
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('ROUGH_WALLS', 'no',  &
                                   val, verbose)
  call String % To_Upper_Case(val)

  select case(val)

    case('YES')
      rough_walls = .true.

    case('NO')
      rough_walls = .false.

    case default
      if(this_proc < 2) then
        print *, '# Unknown wall roughness state :', trim(val)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
