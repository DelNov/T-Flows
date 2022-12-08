!==============================================================================!
  subroutine Connect_Real_Cell(Work,                                    &
                               a01, a02, a03, a04, a05, a06, a07, a08,  &
                               a09, a10, a11, a12, a13, a14, a15, a16)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type)                    :: Work
  real, contiguous, pointer           :: a01(:)
  real, contiguous, pointer, optional :: a02(:), a03(:), a04(:),  &
                                         a05(:), a06(:), a07(:),  &
                                         a08(:), a09(:), a10(:),  &
                                         a11(:), a12(:), a13(:),  &
                                         a14(:), a15(:), a16(:)
!==============================================================================!

  Assert(Work % allocated)

  call Profiler % Start('Work_Mod')

  Work % last_r_cell = Work % last_r_cell + 1
  Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
  a01 => Work % r_cell(Work % last_r_cell) % ptr
  a01(:) = 0.0

  if(present(a02)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a02 => Work % r_cell(Work % last_r_cell) % ptr
    a02(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a03)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a03 => Work % r_cell(Work % last_r_cell) % ptr
    a03(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a04)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a04 => Work % r_cell(Work % last_r_cell) % ptr
    a04(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a05)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a05 => Work % r_cell(Work % last_r_cell) % ptr
    a05(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a06)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a06 => Work % r_cell(Work % last_r_cell) % ptr
    a06(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a07)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a07 => Work % r_cell(Work % last_r_cell) % ptr
    a07(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a08)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a08 => Work % r_cell(Work % last_r_cell) % ptr
    a08(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a09)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a09 => Work % r_cell(Work % last_r_cell) % ptr
    a09(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a10)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a10 => Work % r_cell(Work % last_r_cell) % ptr
    a10(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a11)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a11 => Work % r_cell(Work % last_r_cell) % ptr
    a11(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a12)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a12 => Work % r_cell(Work % last_r_cell) % ptr
    a12(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a13)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a13 => Work % r_cell(Work % last_r_cell) % ptr
    a13(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a14)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a14 => Work % r_cell(Work % last_r_cell) % ptr
    a14(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a15)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a15 => Work % r_cell(Work % last_r_cell) % ptr
    a15(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a16)) then
    Work % last_r_cell = Work % last_r_cell + 1
    Work % max_r_cell  = max(Work % max_r_cell, Work % last_r_cell)
    a16 => Work % r_cell(Work % last_r_cell) % ptr
    a16(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  end subroutine
