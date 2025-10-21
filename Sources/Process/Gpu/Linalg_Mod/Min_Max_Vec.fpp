!==============================================================================!
  subroutine Min_Max_Vec(Lin, n, vmin, vmax, a)
!------------------------------------------------------------------------------!
!>  Find minimum and maximum values in a vector.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)  :: Lin   !! parent class
  integer, intent(in) :: n     !! size of vectors
  real                :: vmin  !! minimum value (result)
  real                :: vmax  !! maximum value (result)
  real                :: a(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  vmin = +HUGE
  vmax = -HUGE
  !$tf-acc loop begin
  do i = 1, n
    vmin = min(vmin, a(i))
    vmax = max(vmax, a(i))
  end do
  !$tf-acc loop end

  end subroutine

