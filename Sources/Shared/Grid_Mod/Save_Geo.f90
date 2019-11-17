!==============================================================================!
  subroutine Grid_Mod_Save_Geo(grid,        &
                               sub,         &  ! subdomain
                               nf_sub,      &  ! number of faces in the sub.
                               nbf_sub)        ! number of buffer cells in sub.
!------------------------------------------------------------------------------!
!   Writes file with geomeatrical data: name.geo                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, nf_sub, nbf_sub
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, s, n, c1, c2, var, subo, fu
  character(len=80)    :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .geo file   !
  !                      !
  !----------------------!
  call File_Mod_Set_Name(name_out, processor=sub, extension='.geo')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

  !----------------------!
  !   Node coordinates   !
  !----------------------!
  do var = 1, 3
    do n=1,grid % n_nodes
      if(grid % new_n(n) > 0) then
        if(var .eq. 1) write(fu) grid % xn(n)
        if(var .eq. 2) write(fu) grid % yn(n)
        if(var .eq. 3) write(fu) grid % zn(n)
      end if
    end do
  end do

  !-----------------------------!
  !   Cell center coordinates   !
  !-----------------------------!
  do var = 1, 3
    do c = 1, grid % n_cells
      if(grid % new_c(c) > 0) then
        if(var .eq. 1) write(fu) grid % xc(c)
        if(var .eq. 2) write(fu) grid % yc(c)
        if(var .eq. 3) write(fu) grid % zc(c)
      end if
    end do
    do s = 1, nbf_sub
      if(var .eq. 1) write(fu) grid % xc(buf_recv_ind(s))
      if(var .eq. 2) write(fu) grid % yc(buf_recv_ind(s))
      if(var .eq. 3) write(fu) grid % zc(buf_recv_ind(s))
    end do
  end do

  !---------------------------!
  !   Boundary cell centers   !
  !---------------------------!
  do var = 1, 3
    do c = -1, -grid % n_bnd_cells, -1
      if(grid % new_c(c) .ne. 0) then
        if(var .eq. 1) write(fu) grid % xc(c)
        if(var .eq. 2) write(fu) grid % yc(c)
        if(var .eq. 3) write(fu) grid % zc(c)
      end if
    end do 
  end do

  !------------------!
  !   Cell volumes   !
  !------------------!
  do c = 1, grid % n_cells
    if(grid % new_c(c) > 0) then
      write(fu) grid % vol(c)
    end if
  end do
  do s = 1, nbf_sub
    write(fu) grid % vol(buf_recv_ind(s))
  end do

  !-------------------!
  !   Wall distance   !
  !-------------------!
  do c = 1, grid % n_cells
    if(grid % new_c(c) > 0) then
      write(fu) grid % wall_dist(c)
    end if
  end do
  do s = 1, nbf_sub
    write(fu) grid % wall_dist(buf_recv_ind(s))
  end do
  do c = -1, -grid % n_bnd_cells, -1
    if(grid % new_c(c) .ne. 0) then
      write(fu) grid % wall_dist(c)
    end if
  end do

  !-----------!
  !   Faces   !
  !-----------!

  ! From 1 to nf_sub -> cell faces for which both cells are inside sub
  do var = 1, 10

    do s = 1, grid % n_faces
      if(grid % new_f(s) > 0 .and. grid % new_f(s) <= nf_sub) then
        if(var .eq.  1)  write(fu) grid % sx(s)
        if(var .eq.  2)  write(fu) grid % sy(s)
        if(var .eq.  3)  write(fu) grid % sz(s)
        if(var .eq.  4)  write(fu) grid % dx(s)
        if(var .eq.  5)  write(fu) grid % dy(s)
        if(var .eq.  6)  write(fu) grid % dz(s)
        if(var .eq.  7)  write(fu) grid % f(s)
        if(var .eq.  8)  write(fu) grid % xf(s)
        if(var .eq.  9)  write(fu) grid % yf(s)
        if(var .eq. 10)  write(fu) grid % zf(s)
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
              if(var .eq.  1)  write(fu) grid % sx(s)
              if(var .eq.  2)  write(fu) grid % sy(s)
              if(var .eq.  3)  write(fu) grid % sz(s)
              if(var .eq.  4)  write(fu) grid % dx(s)
              if(var .eq.  5)  write(fu) grid % dy(s)
              if(var .eq.  6)  write(fu) grid % dz(s)
              if(var .eq.  7)  write(fu) grid % f(s)
              if(var .eq.  8)  write(fu) grid % xf(s)
              if(var .eq.  9)  write(fu) grid % yf(s)
              if(var .eq. 10)  write(fu) grid % zf(s)
            end if  
            if( (grid % comm % proces(c2) .eq. sub) .and.   &
                (grid % comm % proces(c1) .eq. subo) ) then
              if(var .eq.  1)  write(fu) -grid % sx(s)
              if(var .eq.  2)  write(fu) -grid % sy(s)
              if(var .eq.  3)  write(fu) -grid % sz(s)
              if(var .eq.  4)  write(fu) -grid % dx(s)
              if(var .eq.  5)  write(fu) -grid % dy(s)
              if(var .eq.  6)  write(fu) -grid % dz(s)
              if(var .eq.  7)  write(fu) 1.0 - grid % f(s)
              if(var .eq.  8)  write(fu) grid % xf(s) - grid % dx(s)
              if(var .eq.  9)  write(fu) grid % yf(s) - grid % dy(s)
              if(var .eq. 10)  write(fu) grid % zf(s) - grid % dz(s)
            end if
          end if  ! c2 > 0
        end if    ! I think this is not really necessary
      end do  ! s
    end do  ! subo

  end do

  !-----------------!
  !   Periodicity   !
  !-----------------!
  write(fu) grid % per_x
  write(fu) grid % per_y
  write(fu) grid % per_z

  close(fu)

  end subroutine
