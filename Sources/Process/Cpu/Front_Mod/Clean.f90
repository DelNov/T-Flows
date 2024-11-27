!==============================================================================!
  subroutine Clean_Front(Front)
!------------------------------------------------------------------------------!
!>  This subroutine clears the data structures associated with a front.
!>  It resets the front object, removing all elements, vertices, and sides,
!>  preparing it for reuse or disposal.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Deallocates the Elem array, which holds the front elements, and resets   !
!     the number of elements (n_elems) to zero.                                !
!   * Deallocates the Vert array, responsible for storing the vertices of      !
!     the front, and sets the number of vertices (n_verts) to zero.            !
!   * Deallocates the side array, which contains the sides of the front,       !
!     and resets the number of sides (n_sides) to zero.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type),  target :: Front  !! parent class
!==============================================================================!

  !-----------------------------!
  !   Deallocate all elements   !
  !-----------------------------!
  deallocate(Front % Elem)
  Front % n_elems = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(Front % Vert)
  Front % n_verts = 0

  !-----------------------------!
  !   Deallocate all vertices   !
  !-----------------------------!
  deallocate(Front % side)
  Front % n_sides = 0

  end subroutine
