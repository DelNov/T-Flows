!==============================================================================!
  subroutine Surf_Mod_Relax(surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: s, t, e, v, i, j, k
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ns   => surf % n_sides
  ne   => surf % n_elems
  vert => surf % vert
  side => surf % side
  elem => surf % elem

  vert(1:nv) % nne = 0
  do e = 1, ne
    i = elem(e) % i
    j = elem(e) % j
    k = elem(e) % k
    vert(i) % nne = vert(i) % nne + 1
    vert(j) % nne = vert(j) % nne + 1
    vert(k) % nne = vert(k) % nne + 1
  end do

  do t = 6, 3, -1
    do s = 1, ns

      e = vert(side(s) % c) % nne + vert(side(s) % d) % nne  &
        - vert(side(s) % a) % nne - vert(side(s) % b) % nne

      if(e .eq. t) then
        vert(side(s) % a) % nne = vert(side(s) % a) % nne + 1
        vert(side(s) % b) % nne = vert(side(s) % b) % nne + 1
        vert(side(s) % c) % nne = vert(side(s) % c) % nne - 1
        vert(side(s) % d) % nne = vert(side(s) % d) % nne - 1
        call Surf_Mod_Swap_Side(surf, s)
      end if

    end do
  end do

  end subroutine
