!==============================================================================!
  subroutine Find_Front_Elements(Front, verbose)
!------------------------------------------------------------------------------!
!   Finds connectivity for sides and elements                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: cnt_one, cnt_two
  integer                  :: e, eb, ea, c, d, c1, c2, d1, d2, s, n_side
  integer                  :: ss, sum_ijk, sum_cd, i_ver, j_ver, k_ver
  integer                  :: i_s, v1, v2, v3
  integer                  :: ne_tot, ns_tot, cnt_one_tot, cnt_two_tot
  integer, allocatable     :: ci(:), di(:), ei(:), ni(:)
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ns   => Front % n_sides
  ne   => Front % n_elems
  Vert => Front % Vert
  side => Front % side
  Elem => Front % Elem

  !-------------------------------!
  !   Find elements' neighbours   !
  !-------------------------------!
  do s = 1, ns
    c  = side(s) % c
    d  = side(s) % d
    ea = side(s) % ea
    eb = side(s) % eb

    ! Element a
    if(ea > 0) then

      do i_ver = 1, Elem(ea) % nv

        ! Get first and second vertex
        v1 = Elem(ea) % v(i_ver)
        if(i_ver < Elem(ea) % nv) then
          v2 = Elem(ea) % v(i_ver+1)
        else
          v2 = Elem(ea) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          Elem(ea) % ns  = Elem(ea) % ns  + 1
          Elem(ea) % s(Elem(ea) % ns)  = s
          if(eb .gt. 0) then
            Elem(ea) % nne = Elem(ea) % nne + 1
            Elem(ea) % e(Elem(ea) % nne) = eb
          end if
        end if
      end do

    end if  ! ea > 0

    ! Element b
    if(eb > 0) then

      do i_ver = 1, Elem(eb) % nv

        ! Get first and second vertex
        v1 = Elem(eb) % v(i_ver)
        if(i_ver < Elem(eb) % nv) then
          v2 = Elem(eb) % v(i_ver+1)
        else
          v2 = Elem(eb) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          Elem(eb) % ns  = Elem(eb) % ns  + 1
          Elem(eb) % s(Elem(eb) % ns)  = s
          if(ea .gt. 0) then
            Elem(eb) % nne = Elem(eb) % nne + 1
            Elem(eb) % e(Elem(eb) % nne) = ea
          end if
        end if
      end do

    end if  ! eb > 0

  end do

  end subroutine
