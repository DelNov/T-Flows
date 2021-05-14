!==============================================================================!
  subroutine Surf_Mod_Find_Sides(surf, verbose)
!------------------------------------------------------------------------------!
!   Compresses sides' list                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  logical                 :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: cnt_one, cnt_two
  integer                  :: e, eb, ea, i, j, k, b, a, c, d, s, n_side
  integer                  :: ss, sum_ijk, sum_cd, c1, c2, d1, d2
  integer, allocatable     :: ci(:), di(:), ei(:), ni(:)
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ns   => surf % n_sides
  ne   => surf % n_elems
  vert => surf % vert
  side => surf % side
  elem => surf % elem

  !----------------------------!
  !   Put all sides together   !
  !----------------------------!
  n_side = 0  ! initialize side coutner
  do e = 1, ne
    n_side = n_side + 1
    side(n_side) % c  = elem(e) % v(1)
    side(n_side) % d  = elem(e) % v(2)
    side(n_side) % ei = e

    n_side = n_side + 1
    side(n_side) % c  = elem(e) % v(2)
    side(n_side) % d  = elem(e) % v(3)
    side(n_side) % ei = e

    n_side = n_side + 1
    side(n_side) % c  = elem(e) % v(3)
    side(n_side) % d  = elem(e) % v(1)
    side(n_side) % ei = e
  end do
  if(verbose) then
    print *, '# Number of elements:        ', ne
    print *, '# Tentative number of sides: ', n_side
  end if

  allocate(ci(n_side))
  allocate(di(n_side))
  allocate(ei(n_side))
  allocate(ni(n_side))
  do s = 1, n_side
    ci(s) = min(side(s) % c, side(s) % d)
    di(s) = max(side(s) % c, side(s) % d)
    ei(s) = side(s) % ei
    ni(s) = s
  end do

  call Sort % Two_Int_Carry_Int(ci, di, ei)

  !------------------------!
  !   Compress the sides   !
  !------------------------!
  ns = n_side
  n_side = 0

  cnt_one = 0
  cnt_two = 0

  s = 0
  do  ! s = 1, ns

    s = s + 1
    if(s > ns) exit

    ! Take c1, c2, d1 and d2, taking care not to go beyond array boundaries
    c1 = ci(s);  c2 = 0
    d1 = di(s);  d2 = 0
    if(s < ns) then
      c2 = ci(s+1)
      d2 = di(s+1)
    end if

    ! Two sides have the same c1 and c2, handle them both
    if(c1 .eq. c2 .and. d1 .eq. d2) then
      n_side = n_side + 1

      side(n_side) % c = c1  ! here c1 == c2
      side(n_side) % d = d1  ! here d1 == d2

      c = side(n_side) % c
      d = side(n_side) % d

      ! Tedious but (hopefully) correct way to find eb and ea
      do ss = s, s+1

        i = elem(ei(ss)) % v(1)
        j = elem(ei(ss)) % v(2)
        k = elem(ei(ss)) % v(3)

        ! Check ea
        if(i.eq.c .and. j.eq.d) then
          side(n_side) % ea = ei(ss)
          side(n_side) % a  = k
        end if

        if(k.eq.c .and. i.eq.d) then
          side(n_side) % ea = ei(ss)
          side(n_side) % a  = j
        end if

        if(j.eq.c .and. k.eq.d) then
          side(n_side) % ea = ei(ss)
          side(n_side) % a  = i
        end if

        ! Check eb
        if(i.eq.c .and. k.eq.d) then
          side(n_side) % eb = ei(ss)
          side(n_side) % b  = j
        end if

        if(k.eq.c .and. j.eq.d) then
          side(n_side) % eb = ei(ss)
          side(n_side) % b  = i
        end if

        if(j.eq.c .and. i.eq.d) then
          side(n_side) % eb = ei(ss)
          side(n_side) % b  = k
        end if

      end do

      ! You handled two sides, skip one
      cnt_two = cnt_two + 1
      s = s + 1

    ! Side is alone, no twin
    else
      n_side = n_side + 1

      side(n_side) % c = c1
      side(n_side) % d = d1

      c = side(n_side) % c
      d = side(n_side) % d

      ! Tedious but (hopefully) correct way to find eb and ea
      ss = s

      i = elem(ei(ss)) % v(1)
      j = elem(ei(ss)) % v(2)
      k = elem(ei(ss)) % v(3)

      ! Check ea
      if(i.eq.c .and. j.eq.d) then
        side(n_side) % ea = ei(ss)
        side(n_side) % a  = k
      end if

      if(k.eq.c .and. i.eq.d) then
        side(n_side) % ea = ei(ss)
        side(n_side) % a  = j
      end if

      if(j.eq.c .and. k.eq.d) then
        side(n_side) % ea = ei(ss)
        side(n_side) % a  = i
      end if

      ! Check eb
      if(i.eq.c .and. k.eq.d) then
        side(n_side) % eb = ei(ss)
        side(n_side) % b  = j
      end if

      if(k.eq.c .and. j.eq.d) then
        side(n_side) % eb = ei(ss)
        side(n_side) % b  = i
      end if

      if(j.eq.c .and. i.eq.d) then
        side(n_side) % eb = ei(ss)
        side(n_side) % b  = k
      end if

      ! You handled only one side
      cnt_one = cnt_one + 1
    end if
  end do
  ns = n_side
  if(verbose) then
    print *, '# Compressed number of sides:       ', ns
    print *, '# Sides surrounded by two elements: ', cnt_two
    print *, '# Sides surrounded by one element:  ', cnt_one
  end if

  !-------------------------------!
  !   Find elements' neighbours   !
  !-------------------------------!
  do s = 1, ns
    c  = side(s) % c
    d  = side(s) % d
    a  = side(s) % a
    b  = side(s) % b
    eb = side(s) % eb
    ea = side(s) % ea

    ! Element a
    if(ea > 0) then
      if(elem(ea) % v(1) .eq. a) then
        elem(ea) % e(1) = eb
        elem(ea) % s(1) = s
      end if
      if(elem(ea) % v(2) .eq. a) then
        elem(ea) % e(2) = eb
        elem(ea) % s(2) = s
      end if
      if(elem(ea) % v(3) .eq. a) then
        elem(ea) % e(3) = eb
        elem(ea) % s(3) = s
      end if
    end if    ! ea > 0

    ! Element b
    if(eb > 0) then
      if(elem(eb) % v(1) .eq. b) then
        elem(eb) % e(1) = ea
        elem(eb) % s(1) = s
      end if
      if(elem(eb) % v(2) .eq. b) then
        elem(eb) % e(2) = ea
        elem(eb) % s(2) = s
      end if
      if(elem(eb) % v(3) .eq. b) then
        elem(eb) % e(3) = ea
        elem(eb) % s(3) = s
      end if
    end if

  end do

  ! Checking
  do e = 1, ne
    sum_ijk = elem(e) % v(1) + elem(e) % v(2) + elem(e) % v(3)
    sum_cd  = side(elem(e) % s(1)) % c + side(elem(e) % s(1)) % d  &
            + side(elem(e) % s(2)) % c + side(elem(e) % s(2)) % d  &
            + side(elem(e) % s(3)) % c + side(elem(e) % s(3)) % d
    if( sum_cd / sum_ijk .ne. 2 ) then
      print *, '# ERROR in forming elements'' neighbours!'
      stop
    end if
  end do

  call Surf_Mod_Count_Elements_Neighbours(surf)

  end subroutine
