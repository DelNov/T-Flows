!==============================================================================!
  subroutine Surf_Mod_Statistics(surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type),  target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  type(Side_Type), pointer :: side(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: item, i, j, k, c, d, e, s, v, si, sj, sk
  integer                  :: nne_s, nne_e
  real, allocatable        :: nne(:)
  real                     :: a(3), b(3), tri_v(3)
  real                     :: max_rat, min_rat, max_l, min_l
  character(len=160) :: line, n_temp
  integer, parameter :: T=33                ! indent
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ns   => surf % n_sides
  ne   => surf % n_elems
  vert => surf % vert
  side => surf % side
  elem => surf % elem

  !--------------------------!
  !   Compute side lengths   !
  !--------------------------!
  do s = 1, ns
    c = side(s) % c
    d = side(s) % d
    side(s) % length = Math_Mod_Distance(                 &
            vert(c) % x_n, vert(c) % y_n, vert(c) % z_n,  &
            vert(d) % x_n, vert(d) % y_n, vert(d) % z_n)
  end do

  !-----------------------------!
  !   Compute elements' areas   !
  !-----------------------------!
  do e = 1, ne
    i = elem(e) % i
    j = elem(e) % j
    k = elem(e) % k
    a(1) = vert(j) % x_n - vert(i) % x_n
    a(2) = vert(j) % y_n - vert(i) % y_n
    a(3) = vert(j) % z_n - vert(i) % z_n
    b(1) = vert(k) % x_n - vert(i) % x_n
    b(2) = vert(k) % y_n - vert(i) % y_n
    b(3) = vert(k) % z_n - vert(i) % z_n
    tri_v = Math_Mod_Cross_Product(a, b)
    elem(e) % area = sqrt(dot_product(tri_v, tri_v)) * 0.5
  end do

  !------------------------------!
  !   Extreme side size ratios   !
  !------------------------------!
  max_rat = -HUGE
  min_rat = +HUGE
  do e = 1, ne
    si = elem(e) % si
    sj = elem(e) % sj
    sk = elem(e) % sk
    max_l = max(side(si) % length, side(sj) % length, side(sk) % length)
    min_l = min(side(si) % length, side(sj) % length, side(sk) % length)
    max_rat = max(max_rat, max_l/min_l)
    min_rat = min(min_rat, max_l/min_l)
  end do

  !--------------------------------!
  !   Count number of neighbours   !
  !--------------------------------!
  call Surf_Mod_Count_Vertex_Elements(surf)
  nne_s = minval(vert(1:nv) % nne)
  nne_e = maxval(vert(1:nv) % nne)
  allocate(nne(nne_s:nne_e)); nne = 0.0
  do v = 1, nv
    nne(vert(v) % nne) = nne(vert(v) % nne) + 1.0
  end do

  if(this_proc < 2) then

    line( 1:160) = ' '
    line( 1+T:63+T) =   &
               '#=============================================================#'
    print *, trim(line)
    line( 1+T:63+T) =   &
               '#                   Surface mesh statistics                   #'
    print *, trim(line)
    line( 1+T:63+T) =   &
               '#-------------------------------------------------------------#'
    print *, trim(line)

    ! Number of elements (1), nodes (2) and sides (3)
    do item = 1, 3
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      if(item.eq.1) line( 5+T: 5+T+20) = 'Number of elements: '
      if(item.eq.1) write(line(32+T:37+T), '(i6)') surf % n_elems
      if(item.eq.2) line( 5+T: 5+T+20) = 'Number of vertices: '
      if(item.eq.2) write(line(32+T:37+T), '(i6)') surf % n_verts
      if(item.eq.3) line( 5+T: 5+T+20) = 'Number of sides:    '
      if(item.eq.3) write(line(32+T:37+T), '(i6)') surf % n_sides
      print *, trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print *, trim(line)

    ! Maximum (1) and minimum (2) element area
    ! Maximum (3) and minimum (4) side length
    ! Maximum (5) and minimum (6) side length ratio
    do item = 1, 6
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      if(item.eq.1) line( 5+T: 5+T+33) = 'Maximum element area:          '
      if(item.eq.1) write(line(36+T:47+T), '(1pe12.5)')  &
                          maxval(elem(1:ne) % area)
      if(item.eq.2) line( 5+T: 5+T+33) = 'Minimum element area:          '
      if(item.eq.2) write(line(36+T:47+T), '(1pe12.5)')  &
                          minval(elem(1:ne) % area)
      if(item.eq.3) line( 5+T: 5+T+33) = 'Maximum side length:           '
      if(item.eq.3) write(line(36+T:47+T), '(1pe12.5)')  &
                          maxval(side(1:ns) % length)
      if(item.eq.4) line( 5+T: 5+T+33) = 'Minimum side length:           '
      if(item.eq.4) write(line(36+T:47+T), '(1pe12.5)')  &
                          minval(side(1:ns) % length)
      if(item.eq.5) line( 5+T: 5+T+33) = 'Maximum side ratio in element: '
      if(item.eq.5) write(line(36+T:47+T), '(1pe12.5)') max_rat
      if(item.eq.6) line( 5+T: 5+T+33) = 'Minimum side ratio in element: '
      if(item.eq.6) write(line(36+T:47+T), '(1pe12.5)') min_rat
      print *, trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print *, trim(line)

    ! Number of neighbours
    do item = nne_s, nne_e
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      line( 5+T: 5+T+43) = 'Percentage of vertices with XX neighbours: '
      write(line(33+T:34+T), '(i2)') item
      write(line(48+T:53+T), '(f6.2)') nne(item) / nv * 100.0
      line(55+T:55+T) = '%'
      print *, trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print *, trim(line)
    print *, ''

  end if

  end subroutine
