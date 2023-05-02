!==============================================================================!
  integer function Solver_Mod_Linear_Solvers_Code(name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name
!==============================================================================!

  select case(name)

    case('NATIVE')
      Solver_Mod_Linear_Solvers_Code = NATIVE

    case('PETSC')
      Solver_Mod_Linear_Solvers_Code = PETSC

    case default
      call Message % Error(80, 'Unknown linear solver option: '//trim(name)//  &
                               '. This error is critical, exiting! '       //  &
                               'Check the file: Documents/all_control'     //  &
                               '_keywords to see which linear solvers are '//  &
                               'available.  Probably just native and petsc.',  &
                               file=__FILE__, line=__LINE__, one_proc=.true.)
  end select

  end function
