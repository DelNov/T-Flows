!==============================================================================!
  integer function Numerics_Mod_Linear_Solvers_Code(name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name
!==============================================================================!

  select case(name)

    case('NATIVE')
      Numerics_Mod_Linear_Solvers_Code = NATIVE

    case('PETSC')
      Numerics_Mod_Linear_Solvers_Code = PETSC

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown linear solvers option: ',  &
                 trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
      stop

  end select

  end function
