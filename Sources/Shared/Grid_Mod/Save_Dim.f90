!==============================================================================!
  subroutine Save_Dim(Grid, sub)
!------------------------------------------------------------------------------!
!   Writes file with grid dimensions (.dim, used to be .geo)                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: sub
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, s, n, c1, c2, var, fu
  character(SL) :: name_out
!==============================================================================!

  call Profiler % Start('Save_Dim')

  !----------------------!
  !   Create .dim file   !
  !----------------------!
  call File % Set_Name(name_out, processor=sub, extension='.dim')
  call File % Open_For_Writing_Binary(name_out, fu)

  !-------------------------!
  !   Save real precision   !
  !-------------------------!
  write(fu) RP

  !----------------------!
  !   Node coordinates   !
  !----------------------!
  do var = 1, 3
    do n = 1, Grid % n_nodes
      if(Grid % new_n(n) > 0) then
        if(var .eq. 1) write(fu) Grid % xn(n)
        if(var .eq. 2) write(fu) Grid % yn(n)
        if(var .eq. 3) write(fu) Grid % zn(n)
      end if
    end do
  end do

  !-----------------------------!
  !   Cell center coordinates   !
  !-----------------------------!
  do var = 1, 3
    do c = -Grid % n_bnd_cells, Grid % n_cells
      if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
        if(var .eq. 1) write(fu) Grid % xc(Grid % old_c(c))
        if(var .eq. 2) write(fu) Grid % yc(Grid % old_c(c))
        if(var .eq. 3) write(fu) Grid % zc(Grid % old_c(c))
      end if
    end do
  end do

  !-------------------!
  !   Wall distance   !
  !-------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) Grid % wall_dist(Grid % old_c(c))
    end if
  end do

  !------------------!
  !   Cell volumes   !
  !------------------!
  do c = 1, Grid % n_cells
    if(Grid % old_c(c) .ne. 0) then
      write(fu) Grid % vol(Grid % old_c(c))
    end if
  end do

  !--------------------------!
  !   Cell inertia tensors   !
  !--------------------------!
  do var = 1, 6
    do c = 1, Grid % n_cells
      if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
        if(var .eq. 1) write(fu) Grid % ixx(Grid % old_c(c))
        if(var .eq. 2) write(fu) Grid % iyy(Grid % old_c(c))
        if(var .eq. 3) write(fu) Grid % izz(Grid % old_c(c))
        if(var .eq. 4) write(fu) Grid % ixy(Grid % old_c(c))
        if(var .eq. 5) write(fu) Grid % ixz(Grid % old_c(c))
        if(var .eq. 6) write(fu) Grid % iyz(Grid % old_c(c))
      end if
    end do
  end do

  !-----------!
  !   Faces   !
  !-----------!
  do var = 1, 13
    do s = 1, Grid % n_faces + Grid % n_shadows
      if(Grid % old_f(s) .ne. 0) then
        c1 = Grid % faces_c(1, Grid % old_f(s))
        c2 = Grid % faces_c(2, Grid % old_f(s))
        if(Grid % new_c(c2) < 0 .or. Grid % new_c(c1) < Grid % new_c(c2)) then
          if(var .eq.  1)  write(fu) Grid % sx(Grid % old_f(s))
          if(var .eq.  2)  write(fu) Grid % sy(Grid % old_f(s))
          if(var .eq.  3)  write(fu) Grid % sz(Grid % old_f(s))
          if(var .eq.  4)  write(fu) Grid % dx(Grid % old_f(s))
          if(var .eq.  5)  write(fu) Grid % dy(Grid % old_f(s))
          if(var .eq.  6)  write(fu) Grid % dz(Grid % old_f(s))
          if(var .eq.  7)  write(fu) Grid % f (Grid % old_f(s))
          if(var .eq.  8)  write(fu) Grid % xf(Grid % old_f(s))
          if(var .eq.  9)  write(fu) Grid % yf(Grid % old_f(s))
          if(var .eq. 10)  write(fu) Grid % zf(Grid % old_f(s))
          if(var .eq. 11)  write(fu) Grid % rx(Grid % old_f(s))
          if(var .eq. 12)  write(fu) Grid % ry(Grid % old_f(s))
          if(var .eq. 13)  write(fu) Grid % rz(Grid % old_f(s))
        else
          if(var .eq.  1)  write(fu) -Grid % sx(Grid % old_f(s))
          if(var .eq.  2)  write(fu) -Grid % sy(Grid % old_f(s))
          if(var .eq.  3)  write(fu) -Grid % sz(Grid % old_f(s))
          if(var .eq.  4)  write(fu) -Grid % dx(Grid % old_f(s))
          if(var .eq.  5)  write(fu) -Grid % dy(Grid % old_f(s))
          if(var .eq.  6)  write(fu) -Grid % dz(Grid % old_f(s))
          if(var .eq.  7)  write(fu) 1.0 - Grid % f (Grid % old_f(s))
          if(var .eq.  8)  write(fu)  Grid % xf(Grid % old_f(s))
          if(var .eq.  9)  write(fu)  Grid % yf(Grid % old_f(s))
          if(var .eq. 10)  write(fu)  Grid % zf(Grid % old_f(s))
          if(var .eq. 11)  write(fu)  Grid % rx(Grid % old_f(s))
          if(var .eq. 12)  write(fu)  Grid % ry(Grid % old_f(s))
          if(var .eq. 13)  write(fu)  Grid % rz(Grid % old_f(s))
        end if
      end if
    end do
  end do

  !-----------------!
  !   Periodicity   !
  !-----------------!
  write(fu) Grid % per_x
  write(fu) Grid % per_y
  write(fu) Grid % per_z

  close(fu)

  call Profiler % Stop('Save_Dim')

  end subroutine
