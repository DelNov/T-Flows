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
  integer                  :: e, i_ver
!==============================================================================!

  ! Take aliases
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  !---------------------------------!
  !   Browse through all elements   !
  !---------------------------------!
  do e = 1, ne

    elem(e) % xe = 0.0
    elem(e) % ye = 0.0
    elem(e) % ze = 0.0

    do i_ver = 1, elem(e) % nv
      elem(e) % xe = elem(e) % xe + vert(elem(e) % v(i_ver)) % x_n
      elem(e) % ye = elem(e) % ye + vert(elem(e) % v(i_ver)) % y_n
      elem(e) % ze = elem(e) % ze + vert(elem(e) % v(i_ver)) % z_n
    end do

    elem(e) % xe = elem(e) % xe / real(elem(e) % nv)
    elem(e) % ye = elem(e) % ye / real(elem(e) % nv)
    elem(e) % ze = elem(e) % ze / real(elem(e) % nv)

  end do

  end subroutine
