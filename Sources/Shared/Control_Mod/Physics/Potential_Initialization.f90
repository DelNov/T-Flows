!==============================================================================!
  subroutine Control_Mod_Potential_Initialization(pot_init, verbose)
!------------------------------------------------------------------------------!
!   Reading potential initialization from the control file.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: pot_init
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('POTENTIAL_INITIALIZATION', 'no',  &
                                   val, verbose)
  call To_Upper_Case(val)

  select case(val)

    case('YES')
      pot_init = .true.

    case('NO')
      pot_init = .false.

    case default
      if(this_proc < 2) then
        print *, '# Unknown state for potential initalization:', trim(val)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
