!==============================================================================!
  subroutine Surf_Mod_Calculate_Element_Centroids(surf)
!------------------------------------------------------------------------------!
!   Calculates element centroids                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: ne
  integer                  :: e
!==============================================================================!

  ! Take aliases
  ne   => surf % n_elems
  Vert => surf % Vert
  elem => surf % elem

  !---------------------------------!
  !   Browse through all elements   !   (in the future, only for this_proc)
  !---------------------------------!
  do e = 1, ne

    elem(e) % xe = ONE_THIRD * (  Vert(elem(e) % v(1)) % x_n   &
                                + Vert(elem(e) % v(2)) % x_n   &
                                + Vert(elem(e) % v(3)) % x_n )

    elem(e) % ye = ONE_THIRD * (  Vert(elem(e) % v(1)) % y_n   &
                                + Vert(elem(e) % v(2)) % y_n   &
                                + Vert(elem(e) % v(3)) % y_n )

    elem(e) % ze = ONE_THIRD * (  Vert(elem(e) % v(1)) % z_n   &
                                + Vert(elem(e) % v(2)) % z_n   &
                                + Vert(elem(e) % v(3)) % z_n )

  end do

  end subroutine
