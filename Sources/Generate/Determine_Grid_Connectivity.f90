!==============================================================================!
  subroutine Determine_Grid_Connectivity(grid, rrun) 
!------------------------------------------------------------------------------!
!   Determines the topology of the cells, faces and boundary cells.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: rrun
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s
  integer :: c1, c2, m, pass
  integer :: lfn(6,4)
!==============================================================================!

  data    lfn / 1, 1, 2, 4, 3, 5,  &
                2, 5, 6, 8, 7, 7,  &
                4, 6, 8, 7, 5, 8,  &
                3, 2, 4, 3, 1, 6  /

  print *, '# Now determining the topology. This may take a while !'

  !------------------------------!
  !                              !
  !   Count the boundary cells   !
  !                              !
  !------------------------------!
  grid % n_bnd_cells = 0
  do c = 1, grid % n_cells
    do m = 1, 24   ! neighbour cells
      if(grid % cells_c(m,c)  < 0) then
        grid % n_bnd_cells   = grid % n_bnd_cells + 1

        ! Remember the boundary color, take positive value for color
        grid % bnd_cond % color(-grid % n_bnd_cells) =  -grid % cells_c(m,c)  

        ! Put new boundary cell into place  
        grid % cells_c(m,c)  = -grid % n_bnd_cells

        ! Material color
        grid % material(-grid % n_bnd_cells) = grid % material(c)

      end if 
    end do
  end do 

  !---------------------------!
  !   Create the array with   ! 
  !   information on faces    !
  !---------------------------!
  grid % n_faces = 0     ! initialize the number of sides
  do pass = 1,3
    if(pass .eq. 2) WallFacFst = grid % n_faces+1
  do c1=1, grid % n_cells
    do m = 1, 24 ! through all the neighbouring cells
      c2=grid % cells_c(m,c1)
      if( (pass .eq. 1) .and.                             &
          (c2>c1)     .and.                               &
          (grid % material(c1) .eq. grid % material(c2))  &
          .or.                                            &
          (pass .eq. 2) .and.                             &
          (c2>c1)     .and.                               &
          (grid % material(c1) .ne. grid % material(c2))    &
          .or.                                            &
          (pass.eq.3) .and.                               &
          (c2<0) ) then
        grid % n_faces=grid % n_faces+1

        ! Which volumes are connected with side grid % n_faces
        grid % faces_c(1, grid % n_faces)=c1
        grid % faces_c(2,grid % n_faces)=c2 

        ! Which is c2 neighbour of c1 and vice versa
        do c = 1,24
          if(grid % cells_c(c,c1) .eq. c2) then 
            face_c_to_c(grid % n_faces,1) = c
          end if

          if(c2 > 0) then
            if(grid % cells_c(c,c2) .eq. c1) then 
              face_c_to_c(grid % n_faces,2) = c
            end if 
          end if
        end do

        ! Nodes of a side grid % n_faces
        if(c2  > 0) then
          if(level(c2)  > level(c1)) then
            grid % faces_n(1,grid % n_faces) =  &
              grid % cells_n( lfn(face_c_to_c(grid % n_faces,2),4), c2 )
            grid % faces_n(2,grid % n_faces) =  &
              grid % cells_n( lfn(face_c_to_c(grid % n_faces,2),3), c2 )
            grid % faces_n(3,grid % n_faces) =  &
              grid % cells_n( lfn(face_c_to_c(grid % n_faces,2),2), c2 )
            grid % faces_n(4,grid % n_faces) =  &
              grid % cells_n( lfn(face_c_to_c(grid % n_faces,2),1), c2 )
          else
            grid % faces_n(1,grid % n_faces) = grid % cells_n( lfn(m,1), c1 )
            grid % faces_n(2,grid % n_faces) = grid % cells_n( lfn(m,2), c1 )
            grid % faces_n(3,grid % n_faces) = grid % cells_n( lfn(m,3), c1 )
            grid % faces_n(4,grid % n_faces) = grid % cells_n( lfn(m,4), c1 )
          end if
        else
          grid % faces_n(1,grid % n_faces) = grid % cells_n( lfn(m,1), c1 )
          grid % faces_n(2,grid % n_faces) = grid % cells_n( lfn(m,2), c1 )
          grid % faces_n(3,grid % n_faces) = grid % cells_n( lfn(m,3), c1 )
          grid % faces_n(4,grid % n_faces) = grid % cells_n( lfn(m,4), c1 )
        end if 

      end if
    end do   ! m
  end do     ! c1
  end do     ! pass
  WallFacLst = grid % n_faces

  print '(a38,i7)', '# Wall and interface faces start at: ', WallFacFst 
  print '(a38,i7)', '# Wall and interface faces end at  : ', WallFacLst

  if(.not. rrun) then
  grid % n_bnd_cells = 0
  do c = 1, grid % n_cells
    do m = 1, 24   ! neighbour cells
      if(grid % cells_c(m,c)  < 0) then
        grid % n_bnd_cells   = grid % n_bnd_cells + 1

        ! Restore the boundary color, take positive value for color
        grid % cells_c(m,c)  = -grid % bnd_cond % color(-grid % n_bnd_cells)
      end if 
    end do
  end do 
  end if

  if(rrun) then

  !-------------------------------------------!
  !   Find the side oposite on the boundary   !
  !-------------------------------------------!
  ! do s = 1, grid % n_faces
  !   c1 = grid % faces_c(1,s)
  !   c2 = grid % faces_c(2,s)
  !   if(c2  < 0) then
  !     do i=1,6
  !       if(grid % cells_c(i,c1) .eq. c2) then
  !         if(i .eq. 1) grid % faces_c(0,s) = grid % cells_c(6,c1)
  !         if(i .eq. 2) grid % faces_c(0,s) = grid % cells_c(4,c1)
  !         if(i .eq. 3) grid % faces_c(0,s) = grid % cells_c(5,c1)
  !         if(i .eq. 4) grid % faces_c(0,s) = grid % cells_c(2,c1)
  !         if(i .eq. 5) grid % faces_c(0,s) = grid % cells_c(3,c1)
  !         if(i .eq. 6) grid % faces_c(0,s) = grid % cells_c(1,c1)
  !       end if
  !     end do
  !   end if
  ! end do 

  !-------------------------------------------------------!
  !   Find the boundary cell on which you should copy     ! 
  !-------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0 .and. grid % bnd_cond % copy_c(c1) .ne. 0) then
      if(grid % bnd_cond % color(c2) .eq. copy_cond(1,0)) then
         grid % bnd_cond % copy_c(c2) = grid % bnd_cond % copy_c(c1)
      end if
    end if
  end do

  end if ! rrun

  !--------------------------------------!
  !   Is there enough allocated memory   !
  !--------------------------------------!
  if( grid % n_faces  > grid % max_n_faces ) then
    print *, '# Error message from Generator'
    print *, '# The number sides is: ',              grid % n_faces
    print *, '# There is space available only for:', grid % max_n_faces
    print *, '# Increase the number of faces in .dom file'
    print *, '# and recompile the code. Good Luck !'
    stop
  end if 

  end subroutine
