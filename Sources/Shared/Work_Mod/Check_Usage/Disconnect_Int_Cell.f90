!==============================================================================!
  subroutine Disconnect_Int_Cell(Work,                                    &
                                 a01, a02, a03, a04, a05, a06, a07, a08,  &
                                 a09, a10, a11, a12, a13, a14, a15, a16)
!------------------------------------------------------------------------------!
!>  Disconnect_Int_Cell reverses the process of Connect_Int_Cell, releasing
!>  the arrays back to the Work_Type after their temporary usage, ensuring they
!>  are available for reuse elsewhere in the program.  The subroutine assumes
!>  the pointers have previously been connected with Connect_Int_Cell.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type)                       :: Work    !! the singleton Work object
  integer, contiguous, pointer           :: a01(:)  !! integer pointer
  integer, contiguous, pointer, optional :: a02(:), a03(:), a04(:),  &
                                            a05(:), a06(:), a07(:),  &
                                            a08(:), a09(:), a10(:),  &
                                            a11(:), a12(:), a13(:),  &
                                            a14(:), a15(:), a16(:)
    !! additional integer pointer
!==============================================================================!

  call Profiler % Start('Work_Mod')

  Work % last_i_cell = Work % last_i_cell - 1
  nullify(a01)

  if(present(a02)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a02)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a03)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a03)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a04)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a04)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a05)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a05)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a06)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a06)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a07)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a07)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a08)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a08)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a09)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a09)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a10)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a10)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a11)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a11)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a12)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a12)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a13)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a13)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a14)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a14)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a15)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a15)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a16)) then
    Work % last_i_cell = Work % last_i_cell - 1
    nullify(a16)
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  end subroutine
