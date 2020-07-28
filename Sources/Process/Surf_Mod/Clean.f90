!==============================================================================!
  subroutine Surf_Mod_Clean(surf)
!------------------------------------------------------------------------------!
!   Cleans the surface                                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type),  target :: surf
!==============================================================================!

  deallocate(surf % cell_has_vertex)

  ! Dellocate memory for working arrays
  ! (Not used yet, will be used in parallel version)
  ! deallocate(i_work)
  ! deallocate(l_work)
  ! deallocate(r_work)

  !-----------------------------!
  !   Deallocate all elements   !
  !-----------------------------!
  deallocate(surf % elem)
  surf % n_elems = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(surf % vert)
  surf % n_verts = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(surf % side)
  surf % n_sides = 0

  end subroutine
