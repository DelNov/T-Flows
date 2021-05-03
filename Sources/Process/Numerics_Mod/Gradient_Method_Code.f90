!==============================================================================!
  integer function Numerics_Mod_Gradient_Method_Code(name)
!------------------------------------------------------------------------------!
!   Decode the string value on advection scheme from control file              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name
!==============================================================================!

  select case(name)

    case('LEAST_SQUARES')
      Numerics_Mod_Gradient_Method_Code = LEAST_SQUARES
    case('GAUSS_THEOREM')
      Numerics_Mod_Gradient_Method_Code = GAUSS_THEOREM

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown gradient computation method: ', trim(name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End
      stop

  end select

  end function
