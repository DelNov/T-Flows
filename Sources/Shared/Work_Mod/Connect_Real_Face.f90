!==============================================================================!
  subroutine Connect_Real_Face(Work,                                    &
                               a01, a02, a03, a04, a05, a06, a07, a08,  &
                               a09, a10, a11, a12, a13, a14, a15, a16)
!------------------------------------------------------------------------------!
!>  Connects the provided real pointers to the allocated real face
!>  arrays, effectively "borrowing" the space for temporary use.  Once the
!>  connected pointer is not needed any more, it must be released with
!>  Disconnect_Real_Face.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type), target            :: Work    !! the singleton Work object
  real, contiguous, pointer           :: a01(:)  !! real pointer
  real, contiguous, pointer, optional :: a02(:), a03(:), a04(:),  &
                                         a05(:), a06(:), a07(:),  &
                                         a08(:), a09(:), a10(:),  &
                                         a11(:), a12(:), a13(:),  &
                                         a14(:), a15(:), a16(:)
    !! additional real pointer
!==============================================================================!

  call Work % Safe_Connect_Real_Face(a01, a02, a03, a04, a05, a06, a07, a08,  &
                                     a09, a10, a11, a12, a13, a14, a15, a16)

  end subroutine
