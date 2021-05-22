!==============================================================================!
  subroutine Clean_Surf(Surf)
!------------------------------------------------------------------------------!
!   Cleans the surface                                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
!==============================================================================!

  deallocate(Surf % cell_has_vertex)

  ! Dellocate memory for working arrays
  ! (Not used yet, will be used in parallel version)
  ! deallocate(i_work)
  ! deallocate(l_work)
  ! deallocate(r_work)

  !-----------------------------!
  !   Deallocate all elements   !
  !-----------------------------!
  deallocate(Surf % elem)
  Surf % n_elems = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(Surf % Vert)
  Surf % n_verts = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(Surf % side)
  Surf % n_sides = 0

  end subroutine
