!==============================================================================!
  function Les(Turb)
!------------------------------------------------------------------------------!
!   Returns true if any of LES models is engaged                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  logical          :: Les
!==============================================================================!

  ! First assume it is false
  Les = .false.

  ! Then check one by one and set to true if necessary
  if(Turb % model .eq. LES_SMAGORINSKY) then
    Les = .true.
    return
  else if(Turb % model .eq. LES_DYNAMIC) then
    Les = .true.
    return
  else if(Turb % model .eq. LES_WALE) then
    Les = .true.
    return
  else if(Turb % model .eq. LES_TVM) then
    Les = .true.
    return
  end if

  end function
