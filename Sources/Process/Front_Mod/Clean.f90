!==============================================================================!
  subroutine Front_Mod_Clean(front)
!------------------------------------------------------------------------------!
!   Cleans the front                                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type),  target :: front
!==============================================================================!

  deallocate(front % cell_has_vertex)

  ! Dellocate memory for working arrays
  ! (Not used yet, will be used in parallel version)
  ! deallocate(i_work)
  ! deallocate(l_work)
  ! deallocate(r_work)

  !-----------------------------!
  !   Deallocate all elements   !
  !-----------------------------!
  deallocate(front % elem)
  front % n_elems = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(front % vert)
  front % n_verts = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(front % side)
  front % n_sides = 0

  end subroutine
