!==============================================================================!
  subroutine Find_Surf_Elements(Surf, verbose)
!------------------------------------------------------------------------------!
!   Compresses sides' list                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  logical                  :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: cnt_one, cnt_two
  integer                  :: e, eb, ea, i, j, k, c, d, s, n_side
  integer                  :: ss, sum_ijk, sum_cd, i_ver, j_ver, k_ver
  integer                  :: i_s, v1, v2, v3, c1, c2, d1, d2
  integer                  :: ne_tot, ns_tot, cnt_one_tot, cnt_two_tot
  integer, allocatable     :: ci(:), di(:), ei(:), ni(:)
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  Vert => Surf % Vert
  side => Surf % side
  elem => Surf % elem

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
      if(elem(ea) % v(2) .eq. c .and. elem(ea) % v(3) .eq. d  .or. &
         elem(ea) % v(3) .eq. c .and. elem(ea) % v(2) .eq. d) then
        elem(ea) % ns = elem(ea) % ns + 1
        elem(ea) % s(1) = s
        if(eb .gt. 0) then
          elem(ea) % nne = elem(ea) % nne + 1
          elem(ea) % e(1) = eb
        end if
      end if

      if(elem(ea) % v(1) .eq. c .and. elem(ea) % v(3) .eq. d  .or. &
         elem(ea) % v(3) .eq. c .and. elem(ea) % v(1) .eq. d) then
        elem(ea) % ns = elem(ea) % ns + 1
        elem(ea) % s(2) = s
        if(eb .gt. 0) then
          elem(ea) % nne = elem(ea) % nne + 1
          elem(ea) % e(2) = eb
        end if
      end if

      if(elem(ea) % v(1) .eq. c .and. elem(ea) % v(2) .eq. d  .or. &
         elem(ea) % v(2) .eq. c .and. elem(ea) % v(1) .eq. d) then
        elem(ea) % ns = elem(ea) % ns + 1
        elem(ea) % s(3) = s
        if(eb .gt. 0) then
          elem(ea) % nne = elem(ea) % nne + 1
          elem(ea) % e(3) = eb
        end if
      end if
    end if  ! ea > 0

    ! Element b
    if(eb > 0) then
      if(elem(eb) % v(2) .eq. c .and. elem(eb) % v(3) .eq. d  .or. &
         elem(eb) % v(3) .eq. c .and. elem(eb) % v(2) .eq. d) then
        elem(eb) % ns = elem(eb) % ns + 1
        elem(eb) % s(1) = s
        if(ea .gt. 0) then
          elem(eb) % nne = elem(eb) % nne + 1
          elem(eb) % e(1) = ea
        end if
      end if

      if(elem(eb) % v(1) .eq. c .and. elem(eb) % v(3) .eq. d  .or. &
         elem(eb) % v(3) .eq. c .and. elem(eb) % v(1) .eq. d) then
        elem(eb) % ns = elem(eb) % ns + 1
        elem(eb) % s(2) = s
        if(ea .gt. 0) then
          elem(eb) % nne = elem(eb) % nne + 1
          elem(eb) % e(2) = ea
        end if
      end if

      if(elem(eb) % v(1) .eq. c .and. elem(eb) % v(2) .eq. d  .or. &
         elem(eb) % v(2) .eq. c .and. elem(eb) % v(1) .eq. d) then
        elem(eb) % ns = elem(eb) % ns + 1
        elem(eb) % s(3) = s
        if(ea .gt. 0) then
          elem(eb) % nne = elem(eb) % nne + 1
          elem(eb) % e(3) = ea
        end if
      end if
    end if

  end do

  ! Checking
  do e = 1, ne
    sum_ijk = 0
    sum_cd  = 0

    do i_ver = 1, elem(e) % nv
      sum_ijk = sum_ijk + elem(e) % v(i_ver)
    end do

    do i_s = 1, elem(e) % ns
      sum_cd = sum_cd + side(elem(e) % s(i_s)) % c  &
                      + side(elem(e) % s(i_s)) % d
    end do

    if( sum_cd / sum_ijk .ne. 2 ) then
      print *, '# ERROR in forming elements'' neighbours!'
      stop
    end if
  end do

  end subroutine
