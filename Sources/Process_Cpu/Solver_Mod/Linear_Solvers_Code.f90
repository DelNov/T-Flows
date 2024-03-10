!==============================================================================!
  integer function Solver_Mod_Linear_Solvers_Code(name)
!------------------------------------------------------------------------------!
!>  This function in Solver_Mod module converts a textual representation of
!>  the linear solvers used into its corresponding numerical code. The function
!>  ensures that the user input from the control file is accurately translated
!>  to a predefined numerical value that Process can act upon.  The codes for
!>  native and PETSc solvers are defined in Solver_Mod itself.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name  !! string specifying the solvers to use.  It can have
                         !! only two values: 'NATIVE' or 'PETSC'
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
