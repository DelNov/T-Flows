!==============================================================================!
  subroutine Print_Front_Statistics(Front)
!------------------------------------------------------------------------------!
!>  This calculates and displays various statistics of a front mesh.  It is
!>  typically involved at the beginning of a time step in Process, just after
!>  the front has been created which, in turn, is invoked after advancing the
!>  VOF function for the new time step.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Calculation of side lengths for each side in the front.                  !
!   * Determination of extreme side size ratios, such as the maximum and       !
!     minimum ratios.                                                          !
!   * Computing the total surface area of the front.                           !
!   * Counting the number of elements, vertices, and sides, including a        !
!     global sum for parallel execution.                                       !
!   * Finding the maximum and minimum element areas and side lengths.          !
!   * Reporting the percentage of vertices with a certain number of neighbors. !
!   * Outputs these statistics, providing insights into the front mesh's       !
!     geometric properties.                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type),  target :: Front  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  type(Side_Type), pointer :: side(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: nv_tot, item, c, d, e, s, v, si, sj, sk
  integer                  :: nne_s, nne_e
  integer                  :: n_elems_tot, n_verts_tot, n_sides_tot
  real, allocatable        :: nne(:)
  real                     :: max_rat, min_rat, max_l, min_l
  real                     :: min_a, max_a, tot_area
  character(len=160)       :: line
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: T=33  ! indent
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ns   => Front % n_sides
  ne   => Front % n_elems
  Vert => Front % Vert
  side => Front % side
  Elem => Front % Elem

  nv_tot = nv
  if(Front % mesh_divided) then
    call Global % Sum_Int(nv_tot)
  end if

  !--------------------------!
  !   Compute side lengths   !
  !--------------------------!
  do s = 1, ns
    c = side(s) % c
    d = side(s) % d
    side(s) % length = Math % Distance(                   &
            Vert(c) % x_n, Vert(c) % y_n, Vert(c) % z_n,  &
            Vert(d) % x_n, Vert(d) % y_n, Vert(d) % z_n)
  end do

  !------------------------------!
  !   Extreme side size ratios   !
  !------------------------------!
  max_rat = -HUGE
  min_rat = +HUGE
  do e = 1, ne
    si = Elem(e) % s(1)
    sj = Elem(e) % s(2)
    sk = Elem(e) % s(3)
    max_l = max(side(si) % length, side(sj) % length, side(sk) % length)
    min_l = min(side(si) % length, side(sj) % length, side(sk) % length)
    max_rat = max(max_rat, max_l/min_l)
    min_rat = min(min_rat, max_l/min_l)
  end do
  call Global % Max_Real(max_rat)
  call Global % Min_Real(min_rat)

  !------------------------!
  !   Total surface area   !
  !------------------------!
  tot_area = 0.0
  do e = 1, ne
    tot_area = tot_area + Elem(e) % area
  end do
  if(Front % mesh_divided) then
    call Global % Sum_Real(tot_area)
  end if

  !--------------------------------!
  !   Count number of neighbours   !
  !--------------------------------!
  nne_s = minval(Vert(1:nv) % nne)
  nne_e = maxval(Vert(1:nv) % nne)
  call Global % Min_Int(nne_s)
  call Global % Max_Int(nne_e)
  allocate(nne(nne_s:nne_e)); nne = 0.0
  do v = 1, nv
    nne(Vert(v) % nne) = nne(Vert(v) % nne) + 1.0
  end do
  if(Front % mesh_divided) then
    do item = nne_s, nne_e
      call Global % Sum_Real(nne(item))
    end do
  end if

  !------------------------------------------!
  !   Work out values for parellel version   !
  !------------------------------------------!
  n_elems_tot = Front % n_elems
  n_verts_tot = Front % n_verts
  n_sides_tot = Front % n_sides
  if(Front % mesh_divided) then
    call Global % Sum_Int(n_elems_tot)
    call Global % Sum_Int(n_verts_tot)
    call Global % Sum_Int(n_sides_tot)
  end if

  max_l = maxval(side(1:ns) % length)
  min_l = minval(side(1:ns) % length)
  call Global % Max_Real(max_l)
  call Global % Min_Real(min_l)

  max_a = maxval(Elem(1:ne) % area)
  min_a = minval(Elem(1:ne) % area)
  call Global % Max_Real(max_a)
  call Global % Min_Real(min_a)

  if(First_Proc()) then

    line( 1:160) = ' '
    line( 1+T:63+T) =   &
               '#=============================================================#'
    print '(a)', trim(line)
    line( 1+T:63+T) =   &
               '#                  Interface mesh statistics                  #'
    print '(a)', trim(line)
    line( 1+T:63+T) =   &
               '#-------------------------------------------------------------#'
    print '(a)', trim(line)

    ! Number of elements (1), vertices (2) and sides (3)
    do item = 1, 3
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      if(item.eq.1) line( 5+T: 5+T+20) = 'Number of elements: '
      if(item.eq.1) write(line(32+T:37+T), '(i6)') n_elems_tot
      if(item.eq.2) line( 5+T: 5+T+20) = 'Number of vertices: '
      if(item.eq.2) write(line(32+T:37+T), '(i6)') n_verts_tot
      if(item.eq.3) line( 5+T: 5+T+20) = 'Number of sides:    '
      if(item.eq.3) write(line(32+T:37+T), '(i6)') n_sides_tot
      print '(a)', trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print '(a)', trim(line)

    ! Maximum (1) and minimum (2) element area
    ! Maximum (3) and minimum (4) side length
    ! Maximum (5) and minimum (6) side length ratio
    do item = 1, 7
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      if(item.eq.1) line( 5+T: 5+T+33) = 'Maximum element area:          '
      if(item.eq.1) write(line(36+T:47+T), '(1pe12.5)') max_a
      if(item.eq.2) line( 5+T: 5+T+33) = 'Minimum element area:          '
      if(item.eq.2) write(line(36+T:47+T), '(1pe12.5)') min_a
      if(item.eq.3) line( 5+T: 5+T+33) = 'Total surface area  :          '
      if(item.eq.3) write(line(36+T:47+T), '(1pe12.5)')  &
                          tot_area
      if(item.eq.4) line( 5+T: 5+T+33) = 'Maximum side length:           '
      if(item.eq.4) write(line(36+T:47+T), '(1pe12.5)') max_l
      if(item.eq.5) line( 5+T: 5+T+33) = 'Minimum side length:           '
      if(item.eq.5) write(line(36+T:47+T), '(1pe12.5)') min_l
      if(item.eq.6) line( 5+T: 5+T+33) = 'Maximum side ratio in element: '
      if(item.eq.6) write(line(36+T:47+T), '(1pe12.5)') max_rat
      if(item.eq.7) line( 5+T: 5+T+33) = 'Minimum side ratio in element: '
      if(item.eq.7) write(line(36+T:47+T), '(1pe12.5)') min_rat
      print '(a)', trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print '(a)', trim(line)

    ! Number of neighbours
    do item = nne_s, nne_e
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      line( 5+T: 5+T+43) = 'Percentage of vertices with XX neighbours: '
      write(line(33+T:34+T), '(i2)') item
      write(line(48+T:53+T), '(f6.2)') nne(item) / nv_tot * 100.0
      line(55+T:55+T) = '%'
      print '(a)', trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print '(a)', trim(line)
    print '(a)', ''

  end if

  end subroutine
