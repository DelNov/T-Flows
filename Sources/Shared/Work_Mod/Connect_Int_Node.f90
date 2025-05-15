!==============================================================================!
  subroutine Connect_Int_Node(Work,                                    &
                              a01, a02, a03, a04, a05, a06, a07, a08,  &
                              a09, a10, a11, a12, a13, a14, a15, a16)
!------------------------------------------------------------------------------!
!>  Connects the provided integer pointers to the allocated integer node
!>  arrays, effectively "borrowing" the space for temporary use.  Once the
!>  connected pointer is not needed any more, it must be released with
!>  Disconnect_Int_Node.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type),    target            :: Work    !! the singleton Work object
  integer, contiguous, pointer           :: a01(:)  !! integer pointer
  integer, contiguous, pointer, optional :: a02(:), a03(:), a04(:),  &
                                            a05(:), a06(:), a07(:),  &
                                            a08(:), a09(:), a10(:),  &
                                            a11(:), a12(:), a13(:),  &
                                            a14(:), a15(:), a16(:)
    !! additional integer pointer
!==============================================================================!

  call Work % Safe_Connect_Int_Node(a01, a02, a03, a04, a05, a06, a07, a08,  &
                                    a09, a10, a11, a12, a13, a14, a15, a16)

  end subroutine
