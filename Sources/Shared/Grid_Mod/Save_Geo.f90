!==============================================================================!
  subroutine Grid_Mod_Save_Geo(grid,        &
                               sub)
!------------------------------------------------------------------------------!
!   Writes file with geomeatrical data: name.geo                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, s, n, c1, c2, var, subo, fu
  character(SL) :: name_out
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
    do c = -grid % n_bnd_cells, grid % n_cells
      if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
        if(var .eq. 1) write(fu) grid % xc(grid % old_c(c))
        if(var .eq. 2) write(fu) grid % yc(grid % old_c(c))
        if(var .eq. 3) write(fu) grid % zc(grid % old_c(c))
      end if
    end do
  end do

  !-------------------!
  !   Wall distance   !
  !-------------------!
  do c = -grid % n_bnd_cells, grid % n_cells
    if(grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      write(fu) grid % wall_dist(grid % old_c(c))
    end if
  end do

  !------------------!
  !   Cell volumes   !
  !------------------!
  do c = 1, grid % n_cells
    if(grid % old_c(c) .ne. 0) then
      write(fu) grid % vol(grid % old_c(c))
    end if
  end do

  !-----------!
  !   Faces   !
  !-----------!
  do var = 1, 13
    do s = 1, grid % n_faces + grid % n_shadows
      if(grid % old_f(s) .ne. 0) then
        c1 = grid % faces_c(1, grid % old_f(s))
        c2 = grid % faces_c(2, grid % old_f(s))
        if(grid % new_c(c2) < 0 .or. grid % new_c(c1) < grid % new_c(c2)) then
          if(var .eq.  1)  write(fu) grid % sx(grid % old_f(s))
          if(var .eq.  2)  write(fu) grid % sy(grid % old_f(s))
          if(var .eq.  3)  write(fu) grid % sz(grid % old_f(s))
          if(var .eq.  4)  write(fu) grid % dx(grid % old_f(s))
          if(var .eq.  5)  write(fu) grid % dy(grid % old_f(s))
          if(var .eq.  6)  write(fu) grid % dz(grid % old_f(s))
          if(var .eq.  7)  write(fu) grid % f (grid % old_f(s))
          if(var .eq.  8)  write(fu) grid % xf(grid % old_f(s))
          if(var .eq.  9)  write(fu) grid % yf(grid % old_f(s))
          if(var .eq. 10)  write(fu) grid % zf(grid % old_f(s))
          if(var .eq. 11)  write(fu) grid % xr(grid % old_f(s))
          if(var .eq. 12)  write(fu) grid % yr(grid % old_f(s))
          if(var .eq. 13)  write(fu) grid % zr(grid % old_f(s))
        else
          if(var .eq.  1)  write(fu) -grid % sx(grid % old_f(s))
          if(var .eq.  2)  write(fu) -grid % sy(grid % old_f(s))
          if(var .eq.  3)  write(fu) -grid % sz(grid % old_f(s))
          if(var .eq.  4)  write(fu) -grid % dx(grid % old_f(s))
          if(var .eq.  5)  write(fu) -grid % dy(grid % old_f(s))
          if(var .eq.  6)  write(fu) -grid % dz(grid % old_f(s))
          if(var .eq.  7)  write(fu) 1.0 - grid % f (grid % old_f(s))
          if(var .eq.  8)  write(fu) grid % xf(grid % old_f(s))
          if(var .eq.  9)  write(fu) grid % yf(grid % old_f(s))
          if(var .eq. 10)  write(fu) grid % zf(grid % old_f(s))
          if(var .eq. 11)  write(fu) grid % xr(grid % old_f(s))
          if(var .eq. 12)  write(fu) grid % yr(grid % old_f(s))
          if(var .eq. 13)  write(fu) grid % zr(grid % old_f(s))
        end if
      end if
    end do
  end do

  !-----------------!
  !   Periodicity   !
  !-----------------!
  write(fu) grid % per_x
  write(fu) grid % per_y
  write(fu) grid % per_z

  close(fu)

  end subroutine
