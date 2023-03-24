!==============================================================================!
  subroutine Exchange_Real_Array(Comm, length, phi, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm
  integer,          intent(in)    :: length
  real,             intent(inout) :: phi(length)
  integer,          intent(in)    :: dest         ! destination processor
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
  Unused(length)
  Unused(phi)
  Unused(dest)
!==============================================================================!

  end subroutine
