!==============================================================================!
  subroutine Front_Mod_Calculate_Element_Centroids(front)
!------------------------------------------------------------------------------!
!   Calculates element centroids                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: ne
  integer                  :: e
!==============================================================================!

  ! Take aliases
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  !---------------------------------!
  !   Browse through all elements   !   (in the future, only for this_proc)
  !---------------------------------!
  do e = 1, ne

    elem(e) % xe = ONE_THIRD * (  vert(elem(e) % v(1)) % x_n   &
                                + vert(elem(e) % v(2)) % x_n   &
                                + vert(elem(e) % v(3)) % x_n )

    elem(e) % ye = ONE_THIRD * (  vert(elem(e) % v(1)) % y_n   &
                                + vert(elem(e) % v(2)) % y_n   &
                                + vert(elem(e) % v(3)) % y_n )

    elem(e) % ze = ONE_THIRD * (  vert(elem(e) % v(1)) % z_n   &
                                + vert(elem(e) % v(2)) % z_n   &
                                + vert(elem(e) % v(3)) % z_n )

  end do

  end subroutine
