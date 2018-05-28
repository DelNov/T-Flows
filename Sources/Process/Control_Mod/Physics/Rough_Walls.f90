!==============================================================================!
  subroutine Control_Mod_Rough_Walls(verbose)
!-----------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod,      only: YES, NO
  use Turbulence_Mod, only: rough_walls
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('ROUGH_WALLS', 'no',  &
                                   val, verbose)
  call To_Upper_Case(val)

  select case(val)

    case('YES') 
      rough_walls = YES           

    case('NO')                    
      rough_walls = NO        

    case default
      print *, '# Unknown wall roughness state :', trim(val)
      print *, '# Exiting!'
      stop 

  end select

  end subroutine
