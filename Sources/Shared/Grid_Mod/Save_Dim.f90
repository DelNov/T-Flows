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
  integer           :: c, s, n, c1, c2, var, fu, i
  character(SL)     :: name_out
  real, allocatable :: buffer(:)
!==============================================================================!

  call Profiler % Start('Save_Dim')

  !--------------------------!
  !   Allocate local array   !
  !--------------------------!
  allocate(buffer(max( Grid % n_nodes * 3,                            &
                      (Grid % n_bnd_cells + Grid % n_cells + 1) * 3,  &
                       Grid % n_faces + Grid % n_shadows)))

  !----------------------!
  !   Create .dim file   !
  !----------------------!
  call File % Set_Name(name_out, processor=sub, extension='.dim')
  call File % Open_For_Writing_Binary(name_out, fu)

  !-------------------------!
  !   Save real precision   !
  !-------------------------!
  write(fu) RP

  !------------------------------!
  !   Save version of the file   !
  !------------------------------!
  write(fu) VERSION_DIM

  !----------------------!
  !   Node coordinates   !
  !----------------------!
  i = 0  ! you allocated enough memory for all three components, check up
  do var = 1, 3
    do n = 1, Grid % n_nodes
      if(Grid % new_n(n) > 0) then
        i=i+1
        if(var .eq. 1) buffer(i) = Grid % xn(n)
        if(var .eq. 2) buffer(i) = Grid % yn(n)
        if(var .eq. 3) buffer(i) = Grid % zn(n)
      end if
    end do
  end do
  write(fu) buffer(1:i)

  !-----------------------------!
  !   Cell center coordinates   !
  !-----------------------------!
  i = 0  ! you allocated enough memory for all three components, check up
  do var = 1, 3
    do c = -Grid % n_bnd_cells, Grid % n_cells
      if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
        i=i+1
        if(var .eq. 1) buffer(i) = Grid % xc(Grid % old_c(c))
        if(var .eq. 2) buffer(i) = Grid % yc(Grid % old_c(c))
        if(var .eq. 3) buffer(i) = Grid % zc(Grid % old_c(c))
      end if
    end do
  end do
  write(fu) buffer(1:i)

  !-------------------!
  !   Wall distance   !  (buffered together with the next)
  !-------------------!
  i = 0  ! you allocated enough memory for all three components, check up
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      i=i+1;  buffer(i) = Grid % wall_dist(Grid % old_c(c))
    end if
  end do

  !------------------!
  !   Cell volumes   !  (buffered together with the previous)
  !------------------!
  do c = 1, Grid % n_cells
    if(Grid % old_c(c) .ne. 0) then
      i=i+1;  buffer(i) = Grid % vol(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  !--------------------------!
  !   Cell inertia tensors   !
  !--------------------------!
  do var = 1, 6
    i = 0
    do c = 1, Grid % n_cells
      if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
        i=i+1
        if(var .eq. 1) buffer(i) = Grid % ixx(Grid % old_c(c))
        if(var .eq. 2) buffer(i) = Grid % iyy(Grid % old_c(c))
        if(var .eq. 3) buffer(i) = Grid % izz(Grid % old_c(c))
        if(var .eq. 4) buffer(i) = Grid % ixy(Grid % old_c(c))
        if(var .eq. 5) buffer(i) = Grid % ixz(Grid % old_c(c))
        if(var .eq. 6) buffer(i) = Grid % iyz(Grid % old_c(c))
      end if
    end do
    write(fu) buffer(1:i)
  end do

  !-----------!
  !   Faces   !
  !-----------!
  do var = 1, 13
    i = 0
    do s = 1, Grid % n_faces + Grid % n_shadows
      if(Grid % old_f(s) .ne. 0) then
        c1 = Grid % faces_c(1, Grid % old_f(s))
        c2 = Grid % faces_c(2, Grid % old_f(s))
        if(Grid % new_c(c2) < 0 .or. Grid % new_c(c1) < Grid % new_c(c2)) then
          i=i+1
          if(var .eq.  1)  buffer(i) = Grid % sx(Grid % old_f(s))
          if(var .eq.  2)  buffer(i) = Grid % sy(Grid % old_f(s))
          if(var .eq.  3)  buffer(i) = Grid % sz(Grid % old_f(s))
          if(var .eq.  4)  buffer(i) = Grid % dx(Grid % old_f(s))
          if(var .eq.  5)  buffer(i) = Grid % dy(Grid % old_f(s))
          if(var .eq.  6)  buffer(i) = Grid % dz(Grid % old_f(s))
          if(var .eq.  7)  buffer(i) = Grid % f (Grid % old_f(s))
          if(var .eq.  8)  buffer(i) = Grid % xf(Grid % old_f(s))
          if(var .eq.  9)  buffer(i) = Grid % yf(Grid % old_f(s))
          if(var .eq. 10)  buffer(i) = Grid % zf(Grid % old_f(s))
          if(var .eq. 11)  buffer(i) = Grid % rx(Grid % old_f(s))
          if(var .eq. 12)  buffer(i) = Grid % ry(Grid % old_f(s))
          if(var .eq. 13)  buffer(i) = Grid % rz(Grid % old_f(s))
        else
          i=i+1
          if(var .eq.  1)  buffer(i) = -Grid % sx(Grid % old_f(s))
          if(var .eq.  2)  buffer(i) = -Grid % sy(Grid % old_f(s))
          if(var .eq.  3)  buffer(i) = -Grid % sz(Grid % old_f(s))
          if(var .eq.  4)  buffer(i) = -Grid % dx(Grid % old_f(s))
          if(var .eq.  5)  buffer(i) = -Grid % dy(Grid % old_f(s))
          if(var .eq.  6)  buffer(i) = -Grid % dz(Grid % old_f(s))
          if(var .eq.  7)  buffer(i) = 1.0 - Grid % f (Grid % old_f(s))
          if(var .eq.  8)  buffer(i) =  Grid % xf(Grid % old_f(s))
          if(var .eq.  9)  buffer(i) =  Grid % yf(Grid % old_f(s))
          if(var .eq. 10)  buffer(i) =  Grid % zf(Grid % old_f(s))
          if(var .eq. 11)  buffer(i) =  Grid % rx(Grid % old_f(s))
          if(var .eq. 12)  buffer(i) =  Grid % ry(Grid % old_f(s))
          if(var .eq. 13)  buffer(i) =  Grid % rz(Grid % old_f(s))
        end if
      end if
    end do
    write(fu) buffer(1:i)
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
