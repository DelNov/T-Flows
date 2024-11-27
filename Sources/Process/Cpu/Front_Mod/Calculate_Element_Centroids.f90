!==============================================================================!
  subroutine Calculate_Element_Centroids(Front)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for calculating the centroids of elements
!>  in a front.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Iterating through all elements in the front.                             !
!   * Initializing the centroid coordinates (x, y, z) of each element to zero. !
!   * Accumulating the coordinates of each vertex that forms an element.       !
!   * Dividing the total coordinates by the number of vertices to find the     !
!     average, thus determining the centroid's position.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: ne
  integer                  :: e, i_ver
!==============================================================================!

  ! Take aliases
  ne   => Front % n_elems
  Vert => Front % Vert
  Elem => Front % Elem

  !---------------------------------!
  !   Browse through all elements   !
  !---------------------------------!
  do e = 1, ne

    Elem(e) % xe = 0.0
    Elem(e) % ye = 0.0
    Elem(e) % ze = 0.0

    do i_ver = 1, Elem(e) % nv
      Elem(e) % xe = Elem(e) % xe + Vert(Elem(e) % v(i_ver)) % x_n
      Elem(e) % ye = Elem(e) % ye + Vert(Elem(e) % v(i_ver)) % y_n
      Elem(e) % ze = Elem(e) % ze + Vert(Elem(e) % v(i_ver)) % z_n
    end do

    Elem(e) % xe = Elem(e) % xe / real(Elem(e) % nv)
    Elem(e) % ye = Elem(e) % ye / real(Elem(e) % nv)
    Elem(e) % ze = Elem(e) % ze / real(Elem(e) % nv)

  end do

  end subroutine
