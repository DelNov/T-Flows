!==============================================================================!
  subroutine Create_Buffers_And_Save(grid)
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Gen_Mod 
  use Div_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: b, c, n, s, c1, c2, sub, subo, ln
  integer              :: n_nodes_sub, n_cells_sub, n_faces_sub,  &
                          n_bnd_cells_sub, n_buff_sub, NCSsub, n_copy_sub
  character(len=80)    :: name_buf
  integer, allocatable :: side_cell(:,:)
  logical, parameter   :: verbose = .false.
!==============================================================================!
!   Each subdomain needs two buffers: a send buffer and a receive buffer.      !
!   A receive buffer will be stored as aditional boundary cells for each       !
!   subdomain. So each subdomain will have grid % n_bnd_cells physical         !
!   boundary faces and nbb_C-grid % n_bnd_cells buffer bounndary cells.        !
!   It is handy to do it that way, because most of the algorythms can remain   !
!   the same as they are now.  They won't even "know" that they use values     !
!   from other processors.  On the other hand, a sending buffer has to be      !
!   allocated in a new separate array called simply buffer(). An additional    !
!   array is needed to keep track of all the indexes. That one is called       !
!   buffind().  buffind() has stored cell numbers from it's own subdomain      !
!   so that later they can be copied with (well, something like that):         !
!                                                                              !
!   do i=1,BUFFSIZ                                                             !
!     buffer(i) = U(buffind(i))                                                !
!   end do                                                                     !
!------------------------------------------------------------------------------!

  allocate (grid % comm % nbb_s(0 : maxval(grid % comm % proces(:))))
  allocate (grid % comm % nbb_e(0 : maxval(grid % comm % proces(:))))

  !-------------------------------!
  !                               !
  !   Browse through subdomains   !
  !                               !
  !-------------------------------!
  do sub = 1, maxval(grid % comm % proces(:))

    call Name_File(sub, name_buf, '.buf')
    open(9, file=name_buf)
    print *, '# Creating files: ', trim(name_buf)

    write(9,'(a22)') '#--------------------#'
    write(9,'(a22)') '#   Buffer indexes   #'
    write(9,'(a22)') '#--------------------#'

    ! Cells
    n_cells_sub = 0     ! number of cells in subdomain
    do c = 1, grid % n_cells
      new_c(c) = 0
    end do
    do c = 1, grid % n_cells
      if(grid % comm % proces(c) .eq. sub) then
        n_cells_sub = n_cells_sub + 1     ! increase the number of cells in sub.
        new_c(c)    = n_cells_sub         ! assign new (local) cell number 
      end if
    end do

    ! Nodes
    n_nodes_sub = 0     ! number of cells in subdomain
    do n = 1, grid % n_nodes
      new_n(n) = 0
    end do
    do c = 1, grid % n_cells
      if(grid % comm % proces(c) .eq. sub) then
        do ln = 1, grid % cells_n_nodes(c)
          new_n(grid % cells_n(ln,c)) = -1
        end do
      end if
    end do
    do n = 1, grid % n_nodes
      if(new_n(n) .eq. -1) then
        n_nodes_sub = n_nodes_sub+1
        new_n(n) = n_nodes_sub
      end if
    end do

    ! Faces & real boundary cells
    n_faces_sub     = 0  ! number of sides in subdomain
    n_bnd_cells_sub = 0  ! number of real boundary cells in subdomain
    NCSsub = 0
    do s = 1, grid % n_faces
      new_f(s) = 0
    end do
    do c = -grid % n_bnd_cells, -1
      new_c(c) = 0
    end do

    ! Faces step 1/2: inside the domain
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)  
      c2 = grid % faces_c(2,s) 
      if(c2 > 0) then
        if( (grid % comm % proces(c1) .eq. sub) .and.  &
            (grid % comm % proces(c2) .eq. sub) ) then
          n_faces_sub = n_faces_sub+1
          new_f(s) = n_faces_sub
        end if
      end if 
    end do

    ! Faces step 2/2: on the boundaries + bundary cells
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)  
      c2 = grid % faces_c(2,s) 
      if(c2 < 0) then
        if( grid % comm % proces(c1) .eq. sub )  then
          n_faces_sub = n_faces_sub+1
          new_f(s) = n_faces_sub

          n_bnd_cells_sub =  n_bnd_cells_sub + 1  ! increase n. of bnd. cells
          new_c(c2)       = -n_bnd_cells_sub      ! new loc. number of bnd. cell
        end if
      end if 
    end do

    do s = 1, grid % n_copy
      c1 = grid % bnd_cond % copy_s(1,s)
      c2 = grid % bnd_cond % copy_s(2,s)
      if( (grid % comm % proces(c1) .eq. sub) .and.  &
          (grid % comm % proces(c2) .eq. sub) ) then
        NCSsub = NCSsub+1
      end if
    end do

    print '(a,i5,a)', ' #============================================='
    print '(a,i5,a)', ' # Saving subdomain ', sub, ' with:'
    print '(a,i9,a)', ' # ', n_cells_sub, ' cells'
    print '(a,i9,a)', ' # ', n_nodes_sub, ' nodes' 
    print '(a,i9,a)', ' # ', n_faces_sub, ' sides' 
    print '(a,i9,a)', ' # ', n_bnd_cells_sub, ' physical boundary cells' 
    print '(a,i5,a)', ' #---------------------------------------------'

    !--------------------!
    !   Create buffers   !
    !--------------------!
    n_buff_sub = 0
    n_copy_sub = 0
    write(9,'(a30)') '# Number of physical boundary cells:'
    write(9,'(i8)')  n_bnd_cells_sub   
    do subo = 1, maxval(grid % comm % proces(:))
      if(subo .ne. sub) then
        grid % comm % nbb_s(subo) = n_buff_sub + 1

        ! Faces inside the domain
        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)  
          c2 = grid % faces_c(2,s) 
          if(c2  > 0) then
            if( (grid % comm % proces(c1) .eq. sub) .and.  &
                (grid % comm % proces(c2) .eq. subo) ) then
              n_buff_sub = n_buff_sub+1
              buf_send_ind(n_buff_sub) = new_c(c1)  ! buffer send index 
              buf_recv_ind(n_buff_sub) = c2         ! important for coordinate
              buf_pos(n_buff_sub) = -n_bnd_cells_sub-n_buff_sub

              new_f(s) = n_faces_sub+n_buff_sub
            end if
            if( (grid % comm % proces(c2) .eq. sub) .and.  &
                (grid % comm % proces(c1) .eq. subo) ) then
              n_buff_sub = n_buff_sub+1
              buf_send_ind(n_buff_sub) = new_c(c2)  ! buffer send index
              buf_recv_ind(n_buff_sub) = c1         ! important for coordinate
              buf_pos(n_buff_sub) = -n_bnd_cells_sub-n_buff_sub

              new_f(s) = n_faces_sub+n_buff_sub
            end if
          end if  ! c2 > 0
        end do    ! through sides

        ! Faces on the "copy" boundary
        do s = 1, grid % n_copy
          c1 = grid % bnd_cond % copy_s(1,s)  
          c2 = grid % bnd_cond % copy_s(2,s) 
          if( (grid % comm % proces(c1) .eq. sub) .and.  &
              (grid % comm % proces(c2) .eq. subo) ) then
            n_buff_sub = n_buff_sub+1
            n_copy_sub = n_copy_sub+1
            buf_send_ind(n_buff_sub) = new_c(c1) ! buffer send index 
            buf_recv_ind(n_buff_sub) = c2 
            buf_pos(n_buff_sub)= -(-n_bnd_cells_sub-n_buff_sub) ! watch the sign
          end if
          if( (grid % comm % proces(c2) .eq. sub) .and.  &
              (grid % comm % proces(c1) .eq. subo) ) then
            n_buff_sub = n_buff_sub+1
            n_copy_sub = n_copy_sub+1
            buf_send_ind(n_buff_sub) = new_c(c2) ! buffer send index
            buf_recv_ind(n_buff_sub) = c1 
            buf_pos(n_buff_sub)= -(-n_bnd_cells_sub-n_buff_sub) ! watch the sign
          end if
        end do    ! through sides
        grid % comm % nbb_e(subo) = n_buff_sub

        ! Write to buffer file
        write(9,'(A33)') '#-------------------------------#' 
        write(9,'(A33)') '#   Conections with subdomain:  #' 
        write(9,'(A33)') '#-------------------------------#' 
        write(9,'(I8)')  subo 
        write(9,'(A30)') '# Number of local connections:'
        write(9,'(I8)')  grid % comm % nbb_e(subo) -  &
                         grid % comm % nbb_s(subo)+1 
        write(9,'(A37)') '# Local number in a buffer and index:'
        do b = grid % comm % nbb_s(subo),  &
               grid % comm % nbb_e(subo)
          write(9,'(2I8)') b - grid % comm % nbb_s(subo) + 1, buf_send_ind(b) 
        end do
      end if 

    end do ! for subo

    call Save_Cns_Geo(grid,              &
                      sub,               &
                      n_nodes_sub,       &
                      n_cells_sub,       &
                      n_faces_sub,       &
                      n_bnd_cells_sub,   &
                      n_buff_sub,        &
                      n_copy_sub)

    call Save_Vtu_Cells(grid,         &
                        sub,          &
                        n_nodes_sub,  &
                        n_cells_sub)

    call Save_Vtu_Links(grid,             &
                        sub,              &
                        n_nodes_sub,      &
                        n_cells_sub,      &
                        n_faces_sub,      &
                        n_bnd_cells_sub,  &
                        n_buff_sub)

    if(verbose) then
      print '(a)',    ' # Test:'
      print '(a,i9)', ' # n_nodes_sub     :', n_nodes_sub
      print '(a,i9)', ' # n_cells_sub     :', n_cells_sub
      print '(a,i9)', ' # n_faces_sub     :', n_faces_sub
      print '(a,i9)', ' # n_bnd_cells_sub :', n_bnd_cells_sub

      print '(a)',    ' #==============================================='
      print '(a,i9)', ' # Subdomain   ', sub
      print '(a,i9)', ' # Buffer size ', n_buff_sub
      do subo = 1, maxval(grid % comm % proces(:))
        if(subo .ne. sub) then
          print '(a,i9,a,3i9)', ' # Connections with ', subo ,' : ',  &
            grid % comm % nbb_e(subo) -                               &
            grid % comm % nbb_s(subo)+1,                              &
            n_bnd_cells_sub + grid % comm % nbb_s(subo),              &
            n_bnd_cells_sub + grid % comm % nbb_e(subo)
        end if 
      end do ! for subo
      print '(a)',    ' #-----------------------------------------------'
    end if

  end do   ! through subdomains

  close(9)

  !--------------------------------------------------!
  !                                                  !
  !   Save the entire domain with renumbered cells   !
  !                                                  !
  !--------------------------------------------------!
  do n = 1, grid % n_nodes
    new_n(n)=n
  end do
  do c = 1, grid % n_cells
    new_c(c) = 0
  end do
  do s = 1, grid % n_faces
    new_f(s) = 0
  end do

  n_cells_sub = 0     ! number of cells renumbered
  do sub = 1, maxval(grid % comm % proces(:))
    do c = 1, grid % n_cells
      if(grid % comm % proces(c) .eq. sub) then
        n_cells_sub = n_cells_sub+1
        new_c(c) = n_cells_sub
      end if
    end do
  end do

  n_faces_sub = 0     ! number of sides renumbered
  do sub = 1, maxval(grid % comm % proces(:))
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(grid % comm % proces(c1) .eq. sub) then
        n_faces_sub = n_faces_sub+1
        new_f(s) = n_faces_sub
      end if
    end do
  end do
  print '(a,2i9)', ' # Number of faces: ', grid % n_faces, n_faces_sub

  ! It is not sorting nodes ... is it good?  I doubt
  call Grid_Mod_Sort_Cells_By_Index(grid, new_c(1), grid % n_cells)
  call Grid_Mod_Sort_Faces_By_Index(grid, new_f(1), grid % n_faces)

  call Sort_Mod_Int_By_Index(grid % comm % proces(1),  new_c(1), grid % n_cells)
  call Sort_Mod_Int_By_Index(grid % material(1),new_c(1), grid % n_cells)

  ! This is important for plotting the grid with EpsPar()
  call Sort_Mod_Real_By_Index(grid % dx(1), new_f(1), grid % n_faces)
  call Sort_Mod_Real_By_Index(grid % dy(1), new_f(1), grid % n_faces)
  call Sort_Mod_Real_By_Index(grid % dz(1), new_f(1), grid % n_faces)

  allocate(side_cell(grid % n_faces, 2))
  do s = 1, grid % n_faces
    side_cell(s,1) = grid % faces_c(1,s)
    side_cell(s,2) = grid % faces_c(2,s)
  end do
  call Sort_Mod_Int_By_Index(side_cell(1,1), new_f(1), grid % n_faces)
  call Sort_Mod_Int_By_Index(side_cell(1,2), new_f(1), grid % n_faces)
  do s = 1, grid % n_faces
    grid % faces_c(1,s) = side_cell(s,1)
    grid % faces_c(2,s) = side_cell(s,2)
  end do
  deallocate(side_cell)

  call Save_Vtu_Cells(grid, 0, grid % n_nodes, grid % n_cells)
  call Save_Vtu_Faces(grid)

  end subroutine
