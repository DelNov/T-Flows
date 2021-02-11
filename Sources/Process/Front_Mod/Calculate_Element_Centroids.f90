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

    elem(e) % xe = ONE_THIRD * (  vert(elem(e) % i) % x_n   &
                                + vert(elem(e) % j) % x_n   &
                                + vert(elem(e) % k) % x_n )

    elem(e) % ye = ONE_THIRD * (  vert(elem(e) % i) % y_n   &
                                + vert(elem(e) % j) % y_n   &
                                + vert(elem(e) % k) % y_n )

    elem(e) % ze = ONE_THIRD * (  vert(elem(e) % i) % z_n   &
                                + vert(elem(e) % j) % z_n   &
                                + vert(elem(e) % k) % z_n )

  end do

  end subroutine
