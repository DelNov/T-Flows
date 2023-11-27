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

  Work % last_r_face = Work % last_r_face + 1
  Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
  a01 => Work % r_face(Work % last_r_face) % ptr
  a01(:) = 0.0

  if(present(a02)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a02 => Work % r_face(Work % last_r_face) % ptr
    a02(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a03)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a03 => Work % r_face(Work % last_r_face) % ptr
    a03(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a04)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a04 => Work % r_face(Work % last_r_face) % ptr
    a04(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a05)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a05 => Work % r_face(Work % last_r_face) % ptr
    a05(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a06)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a06 => Work % r_face(Work % last_r_face) % ptr
    a06(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a07)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a07 => Work % r_face(Work % last_r_face) % ptr
    a07(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a08)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a08 => Work % r_face(Work % last_r_face) % ptr
    a08(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a09)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a09 => Work % r_face(Work % last_r_face) % ptr
    a09(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a10)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a10 => Work % r_face(Work % last_r_face) % ptr
    a10(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a11)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a11 => Work % r_face(Work % last_r_face) % ptr
    a11(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a12)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a12 => Work % r_face(Work % last_r_face) % ptr
    a12(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a13)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a13 => Work % r_face(Work % last_r_face) % ptr
    a13(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a14)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a14 => Work % r_face(Work % last_r_face) % ptr
    a14(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a15)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a15 => Work % r_face(Work % last_r_face) % ptr
    a15(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  if(present(a16)) then
    Work % last_r_face = Work % last_r_face + 1
    Work % max_r_face  = max(Work % max_r_face, Work % last_r_face)
    a16 => Work % r_face(Work % last_r_face) % ptr
    a16(:) = 0.0
  else
    call Profiler % Stop('Work_Mod')
    return
  end if

  end subroutine
