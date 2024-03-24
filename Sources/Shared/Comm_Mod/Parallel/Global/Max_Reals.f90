!==============================================================================!
  subroutine Max_Reals(Global, phi_01, phi_02, phi_03, phi_04,  &
                               phi_05, phi_06, phi_07, phi_08)
!------------------------------------------------------------------------------!
!>  The subroutine Max_Reals in a parallel computing framework calculates
!>  the global max of up to eight real values across all processors. This
!>  subroutine is vital in scenarios where an aggregate sum of distributed
!>  real value is required for decision-making, post-processing or further
!>  computational tasks.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  real,             intent(inout) :: phi_01  !! global max 1
  real, optional,   intent(inout) :: phi_02  !! global max 2
  real, optional,   intent(inout) :: phi_03  !! global max 3
  real, optional,   intent(inout) :: phi_04  !! global max 4
  real, optional,   intent(inout) :: phi_05  !! global max 5
  real, optional,   intent(inout) :: phi_06  !! global max 6
  real, optional,   intent(inout) :: phi_07  !! global max 7
  real, optional,   intent(inout) :: phi_08  !! global max 8
!-----------------------------------[Locals]-----------------------------------!
  real    :: phi(8)  ! watch out: hard-coded size
  integer :: n, error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
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

  call Global % Max_Real_Array(n, phi(1:n))

  phi_01 = phi(1)
  if(present(phi_02)) phi_02 = phi(2)
  if(present(phi_03)) phi_03 = phi(3)
  if(present(phi_03)) phi_04 = phi(4)
  if(present(phi_05)) phi_05 = phi(5)
  if(present(phi_06)) phi_06 = phi(6)
  if(present(phi_07)) phi_07 = phi(7)
  if(present(phi_08)) phi_08 = phi(8)

  end subroutine
