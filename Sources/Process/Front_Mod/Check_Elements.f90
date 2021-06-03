!==============================================================================!
  subroutine Check_Elements(Front)
!------------------------------------------------------------------------------!
!   Finds connectivity for sides and elements                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
!-----------------------------------[Locals]-----------------------------------!
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  type(Vert_Type), pointer :: Vert(:)
  type(Grid_Type), pointer :: Grid
  integer,         pointer :: ne
  integer                  :: e, sum_ijk, sum_cd, i_ver, i_s, v, c, i_nod, n, s
!==============================================================================!

  ! Take aliases
  ne   => Front % n_elems
  side => Front % side
  Elem => Front % Elem
  Vert => Front % Vert
  Grid => Front % pnt_grid

  ! Checking
  do e = 1, ne
    sum_ijk = 0
    sum_cd  = 0

    do i_ver = 1, Elem(e) % nv
      sum_ijk = sum_ijk + Elem(e) % v(i_ver)
    end do

    do i_s = 1, Elem(e) % ns
      sum_cd = sum_cd + side(Elem(e) % s(i_s)) % c  &
                      + side(Elem(e) % s(i_s)) % d
    end do

    if( sum_cd / sum_ijk .ne. 2 ) then
      print *, '# ERROR in forming elements'' neighbours!'

      print *, '# Element: ', e
      print *, '# Element''s node coordinates read:'
      do i_ver = 1, Elem(e) % nv
        v = Elem(e) % v(i_ver)
        print '(i8, 99es15.5)', v, Vert(v) % x_n, Vert(v) % y_n, Vert(v) % z_n
      end do

      print *, '# Element''s sides and their coordinates: '
      do i_s = 1, Elem(e) % ns
        s = Elem(e) % s(i_s)
        print '(i8, 99es15.5)', s, Vert(side(s) % c) % x_n,  &
                                   Vert(side(s) % c) % y_n,  &
                                   Vert(side(s) % c) % z_n,  &
                                   Vert(side(s) % d) % x_n,  &
                                   Vert(side(s) % d) % y_n,  &
                                   Vert(side(s) % d) % z_n
      end do

      c = Elem(e) % cell
      print *, '# At cell: ', c
      print *, '# Cell''s node coordinates read:'
      do i_nod = 1, Grid % cells_n_nodes(c)
        n = Grid % cells_n(i_nod, c)
        print '(i8, 99es15.5)', n, Grid % xn(n), Grid % yn(n), Grid % zn(n)
      end do

      stop
    end if
  end do

  end subroutine
