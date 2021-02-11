!==============================================================================!
  subroutine Front_Mod_Find_Sides(front, verbose)
!------------------------------------------------------------------------------!
!   Compresses sides' list                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
  logical                 :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: cnt_one, cnt_two
  integer                  :: e, eb, ea, i, j, k, b, a, c, d, s, n_side
  integer                  :: ss, sum_ijk, sum_cd
  integer, allocatable     :: ci(:), di(:), ei(:), ni(:)
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ns   => front % n_sides
  ne   => front % n_elems
  vert => front % vert
  side => front % side
  elem => front % elem

  !----------------------------!
  !   Put all sides together   !
  !----------------------------!
  n_side = 0  ! initialize side coutner
  do e = 1, ne
    n_side = n_side + 1
    side(n_side) % c  = elem(e) % i
    side(n_side) % d  = elem(e) % j
    side(n_side) % ei = e

    n_side = n_side + 1
    side(n_side) % c  = elem(e) % j
    side(n_side) % d  = elem(e) % k
    side(n_side) % ei = e

    n_side = n_side + 1
    side(n_side) % c  = elem(e) % k
    side(n_side) % d  = elem(e) % i
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

  call Sort_Mod_2_Int_Carry_Int(ci, di, ei)

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
    if(s > ns-1) exit

    if(ci(s) .eq. ci(s+1) .and.  &
       di(s) .eq. di(s+1)) then
      n_side = n_side + 1

      side(n_side) % c  = ci(s)
      side(n_side) % d  = di(s)

      c = side(n_side) % c
      d = side(n_side) % d

      ! Tedious but (hopefully) correct way to find eb and ea
      do ss = s, s+1

        i = elem(ei(ss)) % i
        j = elem(ei(ss)) % j
        k = elem(ei(ss)) % k

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

    else
      n_side = n_side + 1

      side(n_side) % c  = ci(s)
      side(n_side) % d  = di(s)

      c = side(n_side) % c
      d = side(n_side) % d

      ! Tedious but (hopefully) correct way to find eb and ea
      ss = s

      i = elem(ei(ss)) % i
      j = elem(ei(ss)) % j
      k = elem(ei(ss)) % k

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
      if(elem(ea) % i .eq. a) then
        elem(ea) % ei = eb
        elem(ea) % si = s
      end if
      if(elem(ea) % j .eq. a) then
        elem(ea) % ej = eb
        elem(ea) % sj = s
      end if
      if(elem(ea) % k .eq. a) then
        elem(ea) % ek = eb
        elem(ea) % sk = s
      end if
    end if    ! ea > 0

    ! Element b
    if(eb > 0) then
      if(elem(eb) % i .eq. b) then
        elem(eb) % ei = ea
        elem(eb) % si = s
      end if
      if(elem(eb) % j .eq. b) then
        elem(eb) % ej = ea
        elem(eb) % sj = s
      end if
      if(elem(eb) % k .eq. b) then
        elem(eb) % ek = ea
        elem(eb) % sk = s
      end if
    end if

  end do

  ! Checking
  do e = 1, ne
    sum_ijk = elem(e) % i + elem(e) % j + elem(e) % k
    sum_cd  = side(elem(e) % si) % c + side(elem(e) % si) % d  &
            + side(elem(e) % sj) % c + side(elem(e) % sj) % d  &
            + side(elem(e) % sk) % c + side(elem(e) % sk) % d
    if( sum_cd / sum_ijk .ne. 2 ) then
      print *, '# ERROR in forming elements'' neighbours!'
      stop
    end if
  end do

  call Front_Mod_Count_Elements_Neighbours(front)

  end subroutine
