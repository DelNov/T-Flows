!==============================================================================!
  real function Reference_Temperature(Flow, t_ref_function)
!------------------------------------------------------------------------------!
!   Imposes variable reference temperature if it is defined by user.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  integer                   :: c
  real                      :: t_ref_function 
!==============================================================================!

  Reference_Temperature = Flow % t_ref

  if(floor(t_ref_function) .ne. 1000000) then
    Reference_Temperature = t_ref_function
  end if

  end function
