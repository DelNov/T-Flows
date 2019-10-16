!==============================================================================!
  subroutine Surf_Mod_Count_Elements_Neighbours(surf)
!------------------------------------------------------------------------------!
!   Counts elements neighbours                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: ne
  integer                  :: e
!==============================================================================!

  ! Take aliases
  ne   => surf % n_elems
  elem => surf % elem

  ! Do the actual count
  elem(1:ne) % nne = 0
  do e = 1, ne
    if(elem(e) % ei .ne. 0) elem(e) % nne = elem(e) % nne + 1
    if(elem(e) % ej .ne. 0) elem(e) % nne = elem(e) % nne + 1
    if(elem(e) % ek .ne. 0) elem(e) % nne = elem(e) % nne + 1
  end do

  end subroutine
