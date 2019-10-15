!==============================================================================!
  subroutine Surf_Mod_Find_Sides(surf)
!------------------------------------------------------------------------------!
!   Compresses nodes' list                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: e, ea, eb, i, j, k, a, b, c, d, s, n_side
  integer                  :: i1, j1, k1, i2, j2, k2, ss, sum_ijk, sum_cd
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
  print *, '# Number of elements:        ', ne
  print *, '# Tentative number of sides: ', n_side

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
  do s = 1, ns, 2

    if(ci(s) .eq. ci(s+1) .and.  &
       di(s) .eq. di(s+1)) then
      n_side = n_side + 1

      side(n_side) % c  = ci(s)
      side(n_side) % d  = di(s)

      c = side(n_side) % c
      d = side(n_side) % d

      ! Tedious but (hopefully) correct way to find ea and eb
      do ss = s, s+1

        i = elem(ei(ss)) % i
        j = elem(ei(ss)) % j
        k = elem(ei(ss)) % k

        ! Check ea
        if(i.eq.c .and. k.eq.d) then
          side(n_side) % ea = ei(ss)
          side(n_side) % a  = j
        end if

        if(k.eq.c .and. j.eq.d) then
          side(n_side) % ea = ei(ss)
          side(n_side) % a  = i
        end if

        if(j.eq.c .and. i.eq.d) then
          side(n_side) % ea = ei(ss)
          side(n_side) % a  = k
        end if

        ! Check eb
        if(i.eq.c .and. j.eq.d) then
          side(n_side) % eb = ei(ss)
          side(n_side) % b  = k
        end if

        if(k.eq.c .and. i.eq.d) then
          side(n_side) % eb = ei(ss)
          side(n_side) % b  = j
        end if

        if(j.eq.c .and. k.eq.d) then
          side(n_side) % eb = ei(ss)
          side(n_side) % b  = i
        end if

      end do

    else
      print *, '# ERROR!  Compression of sides failed!  '
      stop
    end if
  end do
  ns = n_side
  print *, '# Compressed number of sides: ', ns

! ! Checking
! do s = 1, ns
!   WRITE(400, *) '----------------------------------------'
!   WRITE(400, '(5i5)'), s, side(s) % c, side(s) % d, side(s) % a, side(s) % b
!   WRITE(400, '(3i5)'), s, side(s) % ea, side(s) % eb
!   WRITE(400, '(3i5)'), elem(side(s) % ea) % i, elem(side(s) % ea) % j, elem(side(s) % ea) % k
!   WRITE(400, '(3i5)'), elem(side(s) % eb) % i, elem(side(s) % eb) % j, elem(side(s) % eb) % k
! end do

  do s = 1, ns
    c  = side(s) % c
    d  = side(s) % d
    a  = side(s) % a
    b  = side(s) % b
    ea = side(s) % ea
    eb = side(s) % eb

    ! Element a
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

    ! Element b
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
!   WRITE(500, *) '----------------------------------------'
!   WRITE(500, '(3i5)') elem(e) % ei, elem(e) % ej, elem(e) % ek
!   WRITE(500, '(3i5)') elem(e) % i,  elem(e) % j,  elem(e) % k
!   WRITE(500, '(3i5)') side(elem(e) % si) % c, side(elem(e) % si) % d
!   WRITE(500, '(3i5)') side(elem(e) % sj) % c, side(elem(e) % sj) % d
!   WRITE(500, '(3i5)') side(elem(e) % sk) % c, side(elem(e) % sk) % d
!   WRITE(500, *) sum_ijk, sum_cd
  end do

  elem(1:ne) % nne = 0
  do e = 1, ne
    if(elem(e) % ei .ne. 0) elem(e) % nne = elem(e) % nne + 1
    if(elem(e) % ej .ne. 0) elem(e) % nne = elem(e) % nne + 1
    if(elem(e) % ek .ne. 0) elem(e) % nne = elem(e) % nne + 1
  end do

  end subroutine
