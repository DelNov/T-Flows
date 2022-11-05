!==============================================================================!
  subroutine Adjust_Dim(n_new, a)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)                :: n_new        ! new dimension
  integer, allocatable, dimension(:) :: a
!-----------------------------------[Locals]-----------------------------------!
  integer                            :: n_old
  integer, allocatable, dimension(:) :: t
!==============================================================================!

  ! Don't: call Profiler % Start('Adjust_Dim'), it kills performance

  n_old = size(a, 1)
  if(n_new > n_old) then
    allocate(t(n_old))
    t(1:n_old) = a(1:n_old)
    deallocate(a)
    allocate(a(n_new))
    a(1:n_new) = 0
    a(1:n_old) = t(1:n_old)
    deallocate(t)
  end if

  ! Don't: call Profiler % Stop('Adjust_Dim', it kills performance)

  end subroutine
