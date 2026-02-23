!==============================================================================!
  subroutine Sum_Ints(Global, phi_01, phi_02, phi_03, phi_04,  &
                              phi_05, phi_06, phi_07, phi_08)
!------------------------------------------------------------------------------!
!>  The subroutine Sum_Ints in a parallel computing framework calculates
!>  the global sum of up to eight integer values across all processors.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type),  intent(in)    :: Global  !! global communication object
  integer,           intent(inout) :: phi_01  !! global sum 1
  integer, optional, intent(inout) :: phi_02  !! global sum 2
  integer, optional, intent(inout) :: phi_03  !! global sum 3
  integer, optional, intent(inout) :: phi_04  !! global sum 4
  integer, optional, intent(inout) :: phi_05  !! global sum 5
  integer, optional, intent(inout) :: phi_06  !! global sum 6
  integer, optional, intent(inout) :: phi_07  !! global sum 7
  integer, optional, intent(inout) :: phi_08  !! global sum 8
!-----------------------------------[Locals]-----------------------------------!
  integer :: phi(8)  ! watch out: hard-coded size
  integer :: n
!==============================================================================!

  n = 1;  phi(n) = phi_01
  if(present(phi_02)) then
    n = 2;  phi(2) = phi_02
  end if
  if(present(phi_03)) then
    n = 3;  phi(3) = phi_03
  end if
  if(present(phi_04)) then
    n = 4;  phi(4) = phi_04
  end if
  if(present(phi_05)) then
    n = 5;  phi(5) = phi_05
  end if
  if(present(phi_06)) then
    n = 6;  phi(6) = phi_06
  end if
  if(present(phi_07)) then
    n = 7;  phi(7) = phi_07
  end if
  if(present(phi_08)) then
    n = 8;  phi(8) = phi_08
  end if

  call Global % Sum_Int_Array(n, phi(1:n))

  phi_01 = phi(1)
  if(present(phi_02)) phi_02 = phi(2)
  if(present(phi_03)) phi_03 = phi(3)
  if(present(phi_03)) phi_04 = phi(4)
  if(present(phi_05)) phi_05 = phi(5)
  if(present(phi_06)) phi_06 = phi(6)
  if(present(phi_07)) phi_07 = phi(7)
  if(present(phi_08)) phi_08 = phi(8)

  end subroutine
