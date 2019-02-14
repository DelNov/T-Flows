!==============================================================================!
  integer function Numerics_Mod_Time_Integration_Scheme_Code(scheme_name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: scheme_name
!==============================================================================!

  select case(scheme_name)

    case('LINEAR')
      Numerics_Mod_Time_Integration_Scheme_Code = LINEAR
    case('PARABOLIC')
      Numerics_Mod_Time_Integration_Scheme_Code = PARABOLIC

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown time-integration scheme: ',  &
                 trim(scheme_name)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end function
