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
  type(Vert_Type), pointer :: vert(:)
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
  vert => Surf % vert
  side => Surf % side
  elem => Surf % elem

  !----------------------------!
  !   Put all sides together   !
  !----------------------------!
  n_side = 0  ! initialize side coutner
  do e = 1, ne
    do i_ver = 1, elem(e) % nv
      j_ver = i_ver + 1
      if(j_ver > elem(e) % nv) j_ver = 1

      ! Get first and second vertex
      v1 = elem(e) % v(i_ver)
      v2 = elem(e) % v(j_ver)

      n_side = n_side + 1
      side(n_side) % c  = v1
      side(n_side) % d  = v2
      side(n_side) % ei = e
    end do
  end do
  if(verbose) then
    ne_tot = ne
    ns_tot = n_side
    call Comm_Mod_Global_Sum_Int(ne_tot)
    call Comm_Mod_Global_Sum_Int(ns_tot)
    if(this_proc < 2) then
      print *, '# Number of elements:        ', ne_tot
      print *, '# Tentative number of sides: ', ns_tot
    end if
  end if

  if(n_side > 0) then

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

          do i_ver = 1, elem(ei(ss)) % nv
            j_ver = i_ver + 1
            if(j_ver > elem(ei(ss)) % nv) j_ver = 1
            k_ver = j_ver + 1
            if(k_ver > elem(ei(ss)) % nv) k_ver = 1

            ! Get first, second and thrid vertex
            v1 = elem(ei(ss)) % v(i_ver)
            v2 = elem(ei(ss)) % v(j_ver)
            v3 = elem(ei(ss)) % v(k_ver)

            ! Check ea
            if(v1 .eq. c .and. v2 .eq. d) then
              side(n_side) % ea = ei(ss)
              side(n_side) % a  = v3
            end if

            ! Check eb
            if(v2 .eq. c .and. v1 .eq. d) then
              side(n_side) % eb = ei(ss)
              side(n_side) % b  = v3
            end if

          end do

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

        do i_ver = 1, elem(ei(ss)) % nv
          j_ver = i_ver + 1
          if(j_ver > elem(ei(ss)) % nv) j_ver = 1
          k_ver = j_ver + 1
          if(k_ver > elem(ei(ss)) % nv) k_ver = 1

          ! Get first and second vertex
          v1 = elem(ei(ss)) % v(i_ver)
          v2 = elem(ei(ss)) % v(j_ver)
          v3 = elem(ei(ss)) % v(k_ver)

          ! Check ea
          if(v1 .eq. c .and. v2 .eq. d) then
            side(n_side) % ea = ei(ss)
            side(n_side) % a  = v3
          end if

          ! Check eb
          if(v2 .eq. c .and. v1 .eq. d) then
            side(n_side) % eb = ei(ss)
            side(n_side) % b  = v3
          end if

        end do

        ! You handled only one side
        cnt_one = cnt_one + 1
      end if
    end do

  else
    n_side  = 0
    cnt_one = 0
    cnt_two = 0
  end if
  ns = n_side
  if(verbose) then
    ns_tot = ns
    cnt_two_tot = cnt_two
    cnt_one_tot = cnt_one
    call Comm_Mod_Global_Sum_Int(ns_tot)
    call Comm_Mod_Global_Sum_Int(cnt_two_tot)
    call Comm_Mod_Global_Sum_Int(cnt_one_tot)

    if(this_proc < 2) then
      print *, '# Compressed number of sides:       ', ns_tot
      print *, '# Sides surrounded by two elements: ', cnt_two_tot
      print *, '# Sides surrounded by one element:  ', cnt_one_tot
    end if

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
