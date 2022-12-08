!==============================================================================!
  subroutine Connect_Int_Cell(Work,                                    &
                              a01, a02, a03, a04, a05, a06, a07, a08,  &
                              a09, a10, a11, a12, a13, a14, a15, a16)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type)                       :: Work
  integer, contiguous, pointer           :: a01(:)
  integer, contiguous, pointer, optional :: a02(:), a03(:), a04(:),  &
                                            a05(:), a06(:), a07(:),  &
                                            a08(:), a09(:), a10(:),  &
                                            a11(:), a12(:), a13(:),  &
                                            a14(:), a15(:), a16(:)
!==============================================================================!

  Assert(Work % allocated)

  call Profiler % Start('Work_Mod')

  Work % last_i_cell = Work % last_i_cell + 1
  Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
  a01 => Work % i_cell(Work % last_i_cell) % ptr
  a01(:) = 0

  if(present(a02)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a02 => Work % i_cell(Work % last_i_cell) % ptr
    a02(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a03)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a03 => Work % i_cell(Work % last_i_cell) % ptr
    a03(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a04)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a04 => Work % i_cell(Work % last_i_cell) % ptr
    a04(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a05)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a05 => Work % i_cell(Work % last_i_cell) % ptr
    a05(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a06)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a06 => Work % i_cell(Work % last_i_cell) % ptr
    a06(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a07)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a07 => Work % i_cell(Work % last_i_cell) % ptr
    a07(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a08)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a08 => Work % i_cell(Work % last_i_cell) % ptr
    a08(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a09)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a09 => Work % i_cell(Work % last_i_cell) % ptr
    a09(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a10)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a10 => Work % i_cell(Work % last_i_cell) % ptr
    a10(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a11)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a11 => Work % i_cell(Work % last_i_cell) % ptr
    a11(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a12)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a12 => Work % i_cell(Work % last_i_cell) % ptr
    a12(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a13)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a13 => Work % i_cell(Work % last_i_cell) % ptr
    a13(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a14)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a14 => Work % i_cell(Work % last_i_cell) % ptr
    a14(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a15)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a15 => Work % i_cell(Work % last_i_cell) % ptr
    a15(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a16)) then
    Work % last_i_cell = Work % last_i_cell + 1
    Work % max_i_cell  = max(Work % max_i_cell, Work % last_i_cell)
    a16 => Work % i_cell(Work % last_i_cell) % ptr
    a16(:) = 0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  end subroutine
