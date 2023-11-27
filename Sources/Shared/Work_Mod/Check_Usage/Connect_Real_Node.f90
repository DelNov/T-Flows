!==============================================================================!
  subroutine Connect_Real_Node(Work,                                    &
                               a01, a02, a03, a04, a05, a06, a07, a08,  &
                               a09, a10, a11, a12, a13, a14, a15, a16)
!------------------------------------------------------------------------------!
!>  Connects the provided real pointers to the allocated real node
!>  arrays, effectively "borrowing" the space for temporary use.  Once the
!>  connected pointer is not needed any more, it must be released with
!>  Disconnect_Real_Node.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type)                    :: Work    !! the singleton Work object
  real, contiguous, pointer           :: a01(:)  !! real pointer
  real, contiguous, pointer, optional :: a02(:), a03(:), a04(:),  &
                                         a05(:), a06(:), a07(:),  &
                                         a08(:), a09(:), a10(:),  &
                                         a11(:), a12(:), a13(:),  &
                                         a14(:), a15(:), a16(:)
    !! additional real pointer
!==============================================================================!

  Assert(Work % allocated)

  call Profiler % Start('Work_Mod')

  Work % last_r_node = Work % last_r_node + 1
  Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
  a01 => Work % r_node(Work % last_r_node) % ptr
  a01(:) = 0.0

  if(present(a02)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a02 => Work % r_node(Work % last_r_node) % ptr
    a02(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a03)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a03 => Work % r_node(Work % last_r_node) % ptr
    a03(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a04)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a04 => Work % r_node(Work % last_r_node) % ptr
    a04(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a05)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a05 => Work % r_node(Work % last_r_node) % ptr
    a05(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a06)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a06 => Work % r_node(Work % last_r_node) % ptr
    a06(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a07)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a07 => Work % r_node(Work % last_r_node) % ptr
    a07(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a08)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a08 => Work % r_node(Work % last_r_node) % ptr
    a08(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a09)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a09 => Work % r_node(Work % last_r_node) % ptr
    a09(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a10)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a10 => Work % r_node(Work % last_r_node) % ptr
    a10(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a11)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a11 => Work % r_node(Work % last_r_node) % ptr
    a11(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a12)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a12 => Work % r_node(Work % last_r_node) % ptr
    a12(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a13)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a13 => Work % r_node(Work % last_r_node) % ptr
    a13(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a14)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a14 => Work % r_node(Work % last_r_node) % ptr
    a14(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a15)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a15 => Work % r_node(Work % last_r_node) % ptr
    a15(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a16)) then
    Work % last_r_node = Work % last_r_node + 1
    Work % max_r_node  = max(Work % max_r_node, Work % last_r_node)
    a16 => Work % r_node(Work % last_r_node) % ptr
    a16(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  end subroutine
