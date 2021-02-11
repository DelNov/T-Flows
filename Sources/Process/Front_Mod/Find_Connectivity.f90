!==============================================================================!
  subroutine Front_Mod_Find_Connectivity(front, verbose)
!------------------------------------------------------------------------------!
!   Finds connectivity for sides and elements                                  !
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
  integer                  :: e, eb, ea, c, d, s, n_side
  integer                  :: ss, sum_ijk, sum_cd, i_v, i_s, v1, v2
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
    do i_v = 1, elem(e) % nv

      v1 = i_v
      v2 = i_v + 1
      if(v2 > elem(e) % nv) v2 = 1

      n_side = n_side + 1
      side(n_side) % c  = elem(e) % v(v1)
      side(n_side) % d  = elem(e) % v(v2)
      side(n_side) % ei = e
    end do
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

      side(n_side) % c = ci(s)
      side(n_side) % d = di(s)

      c = side(n_side) % c
      d = side(n_side) % d

      ! Tedious but (hopefully) correct way to find eb and ea
      do ss = s, s+1

        do i_v = 1, elem(ei(ss)) % nv

          ! Get first and second vertex
          v1 = elem(ei(ss)) % v(i_v)
          if(i_v < elem(ei(ss)) % nv) then
            v2 = elem(ei(ss)) % v(i_v+1)
          else
            v2 = elem(ei(ss)) % v(1)
          end if

          ! Check ea
          if(v1.eq.c .and. v2.eq.d) then
            side(n_side) % ea = ei(ss)
          end if

          ! Check eb
          if(v2.eq.c .and. v1.eq.d) then
            side(n_side) % eb = ei(ss)
          end if

        end do

      end do

      ! You handled two sides, skip one
      cnt_two = cnt_two + 1
      s = s + 1

    else
      n_side = n_side + 1

      side(n_side) % c = ci(s)
      side(n_side) % d = di(s)

      c = side(n_side) % c
      d = side(n_side) % d

      ! Tedious but (hopefully) correct way to find eb and ea
      ss = s

      do i_v = 1, elem(ei(ss)) % nv

        ! Get first and second vertex
        v1 = elem(ei(ss)) % v(i_v)
        if(i_v < elem(ei(ss)) % nv) then
          v2 = elem(ei(ss)) % v(i_v+1)
        else
          v2 = elem(ei(ss)) % v(1)
        end if

        ! Check ea
        if(v1 .eq. c .and. v2 .eq. d) side(n_side) % ea = ei(ss)

        ! Check eb
        if(v2 .eq. c .and. v1 .eq. d) side(n_side) % eb = ei(ss)

      end do

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
    ea = side(s) % ea
    eb = side(s) % eb

    ! Element a
    if(ea > 0) then

      do i_v = 1, elem(ea) % nv

        ! Get first and second vertex
        v1 = elem(ea) % v(i_v)
        if(i_v < elem(ea) % nv) then
          v2 = elem(ea) % v(i_v+1)
        else
          v2 = elem(ea) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          elem(ea) % nne = elem(ea) % nne + 1
          elem(ea) % ns  = elem(ea) % ns  + 1
          elem(ea) % e(elem(ea) % nne) = eb
          elem(ea) % s(elem(ea) % ns)  = s
        end if
      end do

    end if  ! ea > 0

    ! Element b
    if(eb > 0) then

      do i_v = 1, elem(eb) % nv

        ! Get first and second vertex
        v1 = elem(eb) % v(i_v)
        if(i_v < elem(eb) % nv) then
          v2 = elem(eb) % v(i_v+1)
        else
          v2 = elem(eb) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          elem(eb) % nne = elem(eb) % nne + 1
          elem(eb) % ns  = elem(eb) % ns  + 1
          elem(eb) % e(elem(eb) % nne) = ea
          elem(eb) % s(elem(eb) % ns)  = s
        end if
      end do

    end if  ! eb > 0

  end do

  ! Checking
  do e = 1, ne
    sum_ijk = 0
    sum_cd  = 0

    do i_v = 1, elem(e) % nv
      sum_ijk = sum_ijk + elem(e) % v(i_v)
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
