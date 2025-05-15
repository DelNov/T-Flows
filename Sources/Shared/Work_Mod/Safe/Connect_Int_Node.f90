!==============================================================================!
  subroutine Safe_Connect_Int_Node(Work,                                    &
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

  call Profiler % Start('Work_Mod')

  Work % last_i_node = Work % last_i_node + 1
  if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
    call Work % Allocate_Int_Node(Work % last_i_node)
  Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
  a01 => Work % i_node(Work % last_i_node) % array
  a01(:) = 0

  if(present(a02)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a02 => Work % i_node(Work % last_i_node) % array
    a02(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a03)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a03 => Work % i_node(Work % last_i_node) % array
    a03(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a04)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a04 => Work % i_node(Work % last_i_node) % array
    a04(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a05)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a05 => Work % i_node(Work % last_i_node) % array
    a05(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a06)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a06 => Work % i_node(Work % last_i_node) % array
    a06(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a07)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a07 => Work % i_node(Work % last_i_node) % array
    a07(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a08)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a08 => Work % i_node(Work % last_i_node) % array
    a08(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a09)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a09 => Work % i_node(Work % last_i_node) % array
    a09(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a10)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a10 => Work % i_node(Work % last_i_node) % array
    a10(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a11)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a11 => Work % i_node(Work % last_i_node) % array
    a11(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a12)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a12 => Work % i_node(Work % last_i_node) % array
    a12(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a13)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a13 => Work % i_node(Work % last_i_node) % array
    a13(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a14)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a14 => Work % i_node(Work % last_i_node) % array
    a14(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a15)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a15 => Work % i_node(Work % last_i_node) % array
    a15(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a16)) then
    Work % last_i_node = Work % last_i_node + 1
    if(.not. allocated(Work % i_node(Work % last_i_node) % array))  &
      call Work % Allocate_Int_Node(Work % last_i_node)
    Work % max_i_node  = max(Work % max_i_node, Work % last_i_node)
    a16 => Work % i_node(Work % last_i_node) % array
    a16(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  end subroutine
