!==============================================================================!
  subroutine Clean_Surf(Surf)
!------------------------------------------------------------------------------!
!>  This subroutine clears the data structures associated with a surface.
!>  It resets the surface object, removing all elements, vertices, and sides,
!>  preparing it for reuse or disposal.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Deallocates the Elem array, which holds the surface elements, and        !
!     resets the number of elements (n_elems) to zero.                         !
!   * Deallocates the Vert array, responsible for storing the vertices of      !
!     the surface, and sets the number of vertices (n_verts) to zero.          !
!   * Deallocates the side array, which contains the sides of the surface,     !
!     and resets the number of sides (n_sides) to zero.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
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
  deallocate(Surf % Elem)
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
