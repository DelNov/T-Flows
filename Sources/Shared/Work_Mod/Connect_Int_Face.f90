!==============================================================================!
  subroutine Connect_Int_Face(Work,                                    &
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

  Work % last_i_face = Work % last_i_face + 1
  Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
  a01 => Work % i_face(Work % last_i_face) % ptr

  if(present(a02)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a02 => Work % i_face(Work % last_i_face) % ptr
    a02(:) = 0
  else
    return
  end if

  if(present(a03)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a03 => Work % i_face(Work % last_i_face) % ptr
    a03(:) = 0
  else
    return
  end if

  if(present(a04)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a04 => Work % i_face(Work % last_i_face) % ptr
    a04(:) = 0
  else
    return
  end if

  if(present(a05)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a05 => Work % i_face(Work % last_i_face) % ptr
    a05(:) = 0
  else
    return
  end if

  if(present(a06)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a06 => Work % i_face(Work % last_i_face) % ptr
    a06(:) = 0
  else
    return
  end if

  if(present(a07)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a07 => Work % i_face(Work % last_i_face) % ptr
    a07(:) = 0
  else
    return
  end if

  if(present(a08)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a08 => Work % i_face(Work % last_i_face) % ptr
    a08(:) = 0
  else
    return
  end if

  if(present(a09)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a09 => Work % i_face(Work % last_i_face) % ptr
    a09(:) = 0
  else
    return
  end if

  if(present(a10)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a10 => Work % i_face(Work % last_i_face) % ptr
    a10(:) = 0
  else
    return
  end if

  if(present(a11)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a11 => Work % i_face(Work % last_i_face) % ptr
    a11(:) = 0
  else
    return
  end if

  if(present(a12)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a12 => Work % i_face(Work % last_i_face) % ptr
    a12(:) = 0
  else
    return
  end if

  if(present(a13)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a13 => Work % i_face(Work % last_i_face) % ptr
    a13(:) = 0
  else
    return
  end if

  if(present(a14)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a14 => Work % i_face(Work % last_i_face) % ptr
    a14(:) = 0
  else
    return
  end if

  if(present(a15)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a15 => Work % i_face(Work % last_i_face) % ptr
    a15(:) = 0
  else
    return
  end if

  if(present(a16)) then
    Work % last_i_face = Work % last_i_face + 1
    Work % max_i_face  = max(Work % max_i_face, Work % last_i_face)
    a16 => Work % i_face(Work % last_i_face) % ptr
    a16(:) = 0
  else
    return
  end if

  end subroutine
