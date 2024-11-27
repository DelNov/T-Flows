!==============================================================================!
  function Rsm(Turb)
!------------------------------------------------------------------------------!
!   Returns true if any of LES models is engaged                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  logical          :: Rsm
!==============================================================================!

  ! First assume it is false
  Rsm = .false.

  ! Then check one by one and set tu tre if necessary
  if(Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    Rsm = .true.
    return
  else if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
    Rsm = .true.
    return
  end if

  end function
