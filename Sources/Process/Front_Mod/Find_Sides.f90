!==============================================================================!
  subroutine Find_Sides(Front, verbose)
!------------------------------------------------------------------------------!
!   Finds connectivity for sides and elements                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
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
  vert => Front % vert
  side => Front % side
  elem => Front % elem

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

  end subroutine
