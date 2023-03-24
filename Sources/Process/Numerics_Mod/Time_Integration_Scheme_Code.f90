!==============================================================================!
  integer function Numerics_Mod_Time_Integration_Scheme_Code(name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: name
!==============================================================================!

  select case(name)

    case('LINEAR')
      Numerics_Mod_Time_Integration_Scheme_Code = LINEAR
    case('PARABOLIC')
      Numerics_Mod_Time_Integration_Scheme_Code = PARABOLIC
    case('RUNGE_KUTTA_3')
      Numerics_Mod_Time_Integration_Scheme_Code = RUNGE_KUTTA_3

    case default
      call Message % Error(60,                                    &
               'Unknown time-integration scheme: '//trim(name)//  &
               '. \n Exiting!', one_proc=.true.)
  end select

  end function
