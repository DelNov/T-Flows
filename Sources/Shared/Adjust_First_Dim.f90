!==============================================================================!
  subroutine Adjust_First_Dim(n_new, a)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)                  :: n_new        ! new dimension
  integer, allocatable, dimension(:,:) :: a
!-----------------------------------[Locals]-----------------------------------!
  integer                              :: n_old, l_old, u_old
  integer, allocatable, dimension(:,:) :: t
!==============================================================================!

  ! Don't: call Profiler % Start('Adjust_First_Dim'), it kills performance

  n_old = size  (a, 1)
  l_old = lbound(a, 2)
  u_old = ubound(a, 2)
  if(n_new > n_old) then
    allocate(t(n_old, l_old:u_old))
    t(1:n_old, l_old:u_old) = a(1:n_old, l_old:u_old)
    deallocate(a)
    allocate(a(n_new, l_old:u_old))
    a(1:n_new, l_old:u_old) = 0
    a(1:n_old, l_old:u_old) = t(1:n_old, l_old:u_old)
    deallocate(t)
  end if

  ! Don't: call Profiler % Stop('Adjust_First_Dim', it kills performance)

  end subroutine
