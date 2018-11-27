!==============================================================================!
  subroutine Save_Geo(grid,        &
                      sub,         &  ! subdomain
                      nf_sub,      &  ! number of faces in the sub.
                      nbf_sub)        ! number of buffer cells in sub.
!------------------------------------------------------------------------------!
!   Writes: name.geo                                                           !
!----------------------------------[Modules]-----------------------------------!
  use Div_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, nf_sub, nbf_sub
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, s, n, c1, c2, count, var, subo
  integer              :: lower_bound, upper_bound
  character(len=80)    :: name_out
  integer, allocatable :: iwork(:,:)
  real, allocatable    :: work(:)
!==============================================================================!
!   sub     - subdomain number                                                 !
!   nf_sub  - number of faces in subdomain, but without faces on buffer        !
!   nbf_sub - number of buffer boundary faces in subdomain                     !
!------------------------------------------------------------------------------!

  lower_bound = min(-nbf_sub, -grid % n_bnd_cells)
  upper_bound = max(grid % n_cells*8, grid % n_faces*4)

  allocate(iwork(lower_bound:upper_bound, 0:2));  iwork = 0
  allocate(work(grid % n_faces));                  work = 0.

  !----------------------!
  !                      !
  !   Create .geo file   !
  !                      !
  !----------------------!
  call Name_File( sub, name_out, '.geo' )
  open(9, file=name_out, form='unformatted', access='stream')
  write(*, *) '# Creating the file: ', trim(name_out)

  !----------------------!
  !   Node coordinates   !
  !----------------------!
  do var = 1, 3
    count = 0
    do n=1,grid % n_nodes
      if(grid % new_n(n) > 0) then
        count = count + 1
        if(var .eq. 1) work(count) = grid % xn(n)
        if(var .eq. 2) work(count) = grid % yn(n)
        if(var .eq. 3) work(count) = grid % zn(n)
      end if
    end do 
    write(9) (work(n), n=1,count)
  end do

  !-----------------------------!
  !   Cell center coordinates   !
  !-----------------------------!
  do var = 1, 3
    count = 0
    do c = 1, grid % n_cells
      if(grid % new_c(c) > 0) then
        count = count + 1
        if(var .eq. 1) work(count) = grid % xc(c)
        if(var .eq. 2) work(count) = grid % yc(c)
        if(var .eq. 3) work(count) = grid % zc(c)
      end if
    end do
    do s = 1, nbf_sub
      count = count + 1
      if(var .eq.  1) work(count) = grid % xc(buf_recv_ind(s))
      if(var .eq.  2) work(count) = grid % yc(buf_recv_ind(s))
      if(var .eq.  3) work(count) = grid % zc(buf_recv_ind(s))
    end do

    write(9) (work(c), c = 1, count)
  end do

  !---------------------------!
  !   Boundary cell centers   !
  !---------------------------!
  do var = 1, 3
    count = 0
    do c = -1, -grid % n_bnd_cells, -1
      if(grid % new_c(c) .ne. 0) then
        count = count + 1
        if(var .eq. 1) work(count) = grid % xc(c)
        if(var .eq. 2) work(count) = grid % yc(c)
        if(var .eq. 3) work(count) = grid % zc(c)
      end if
    end do 
    write(9) (work(c), c = 1, count)
  end do

  !------------------!
  !   Cell volumes   !
  !------------------!
  count = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) > 0) then
      count = count + 1
      work(count) = grid % vol(c)
    end if
  end do
  do s = 1, nbf_sub
    count = count + 1
    work(count) = grid % vol(buf_recv_ind(s))
  end do
  write(9) (work(c), c = 1, count)

  !----------------!
  !   Cell delta   !
  !----------------!
  count = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) > 0) then
      count = count + 1
      work(count) = grid % delta(c)
    end if
  end do
  do s = 1, nbf_sub
    count = count + 1
    work(count) = grid % delta(buf_recv_ind(s))
  end do
  write(9) (work(c), c = 1, count)

  !-------------------!
  !   Wall distance   !
  !-------------------!
  count = 0
  do c = 1, grid % n_cells
    if(grid % new_c(c) > 0) then
      count = count + 1
      work(count) = grid % wall_dist(c)
    end if
  end do
  do s = 1, nbf_sub
    count = count + 1
    work(count) = grid % wall_dist(buf_recv_ind(s))
  end do

  write(9) (work(c), c = 1, count)

  !-----------!
  !   Faces   !
  !-----------!

  ! From 1 to nf_sub -> cell faces for which both cells are inside sub
  do var = 1, 10
  count = 0

  do s = 1, grid % n_faces
    if(grid % new_f(s) > 0 .and. grid % new_f(s) <= nf_sub) then
      count = count + 1
      if(var .eq.  1)  work(count) = grid % sx(s)
      if(var .eq.  2)  work(count) = grid % sy(s)
      if(var .eq.  3)  work(count) = grid % sz(s)
      if(var .eq.  4)  work(count) = grid % dx(s)
      if(var .eq.  5)  work(count) = grid % dy(s)
      if(var .eq.  6)  work(count) = grid % dz(s)
      if(var .eq.  7)  work(count) = grid % f(s)
      if(var .eq.  8)  work(count) = grid % xf(s)
      if(var .eq.  9)  work(count) = grid % yf(s)
      if(var .eq. 10)  work(count) = grid % zf(s)
    end if
  end do

  ! From nf_sub+1 to nf_sub + nbf_sub
  ! (think: are they in right order ?)
  do subo = 1, maxval(grid % comm % proces(:))
    do s = 1, grid % n_faces
      if(grid % new_f(s) > nf_sub .and.  &
         grid % new_f(s) <= nf_sub + nbf_sub) then
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 > 0) then
          if( (grid % comm % proces(c1) .eq. sub) .and.  &
              (grid % comm % proces(c2) .eq. subo) ) then
            count = count + 1
            if(var .eq.  1)  work(count) = grid % sx(s)
            if(var .eq.  2)  work(count) = grid % sy(s)
            if(var .eq.  3)  work(count) = grid % sz(s)
            if(var .eq.  4)  work(count) = grid % dx(s)
            if(var .eq.  5)  work(count) = grid % dy(s)
            if(var .eq.  6)  work(count) = grid % dz(s)
            if(var .eq.  7)  work(count) = grid % f(s)
            if(var .eq.  8)  work(count) = grid % xf(s)
            if(var .eq.  9)  work(count) = grid % yf(s)
            if(var .eq. 10)  work(count) = grid % zf(s)
          end if  
          if( (grid % comm % proces(c2) .eq. sub) .and.   &
              (grid % comm % proces(c1) .eq. subo) ) then
            count = count + 1
            if(var .eq.  1)  work(count) = -grid % sx(s)
            if(var .eq.  2)  work(count) = -grid % sy(s)
            if(var .eq.  3)  work(count) = -grid % sz(s)
            if(var .eq.  4)  work(count) = -grid % dx(s)
            if(var .eq.  5)  work(count) = -grid % dy(s)
            if(var .eq.  6)  work(count) = -grid % dz(s)
            if(var .eq.  7)  work(count) = 1.0 - grid % f(s)
            if(var .eq.  8)  work(count) = grid % xf(s) - grid % dx(s)
            if(var .eq.  9)  work(count) = grid % yf(s) - grid % dy(s)
            if(var .eq. 10)  work(count) = grid % zf(s) - grid % dz(s)
          end if
        end if  ! c2 > 0
      end if    ! I think this is not really necessary
    end do
  end do

  write(9) (work(s), s = 1, count)

  end do

  close(9)

  deallocate (iwork)
  deallocate (work)

  end subroutine
