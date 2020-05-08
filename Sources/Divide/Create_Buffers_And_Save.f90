!==============================================================================!
  subroutine Create_Buffers_And_Save(grid)
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
  use Div_Mod
  use Grid_Mod,      only: Grid_Type,                     &
                           Grid_Mod_Sort_Cells_By_Index,  &
                           Grid_Mod_Sort_Faces_By_Index,  &
                           Grid_Mod_Save_Cns,             &
                           Grid_Mod_Save_Geo
  use Sort_Mod       ! it's a collection of subroutines, no need for "only"
  use Save_Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: b, c, n, s, c1, c2, sub, subo, ln, fu
  integer              :: nn_sub,    &  ! number of nodes in the subdomain
                          nc_sub,    &  ! number of cells in the subdomain
                          nf_sub,    &  ! number of faces in the subdomain
                          nbc_sub,   &  ! number of boundary cells in sub.
                          nbf_sub,   &
                          nbfc_sub
  character(len=80)    :: name_buf
  integer, allocatable :: side_cell(:,:)
  logical, parameter   :: verbose = .false.
!==============================================================================!
!   Each subdomain needs two buffers: a send buffer and a receive buffer.      !
!   A receive buffer will be stored as aditional cells for each subdomain.     !
!   So each subdomain will have grid % n_cells cells, which entials physical   !
!   and buffer cells. It is handy to do it that way, because most of the       !
!   algorythms can remain the same as in sequential run. On the other hand,    !
!   a sending buffer has to be allocated in a new separate array called        !
!   simply buffer(). An additional array is needed to keep track of all the    !
!   indexes. That one is called buffind().                                     !
!------------------------------------------------------------------------------!

  allocate(grid % comm % buff_s_cell(0 : maxval(grid % comm % cell_proc(:))))
  allocate(grid % comm % buff_e_cell(0 : maxval(grid % comm % cell_proc(:))))

  !-------------------------------!
  !                               !
  !   Browse through subdomains   !
  !                               !
  !-------------------------------!
  do sub = 1, maxval(grid % comm % cell_proc(:))

    if(verbose) then
      call File_Mod_Set_Name(name_buf, processor=sub, extension='.buf')
      call File_Mod_Open_File_For_Writing(name_buf, fu)

      write(fu,'(a22)') '#--------------------#'
      write(fu,'(a22)') '#   Buffer indexes   #'
      write(fu,'(a22)') '#--------------------#'
    else
      print *, '# Creating buffers '
    end if

    !-----------!
    !   Cells   !
    !-----------!
    nc_sub = 0     ! number of cells in subdomain
    do c = 1, grid % n_cells
      grid % new_c(c) = 0
    end do
    do c = 1, grid % n_cells
      if(grid % comm % cell_proc(c) .eq. sub) then
        nc_sub   = nc_sub + 1     ! increase the number of cells in sub.
        grid % new_c(c) = nc_sub         ! assign new (local) cell number 
      end if
    end do

    !-----------!
    !   Faces   !
    !-----------!

    ! Faces & real boundary cells
    nf_sub  = 0  ! number of faces in subdomain
    nbc_sub = 0  ! number of real boundary cells in subdomain
    do s = 1, grid % n_faces
      grid % new_f(s) = 0
    end do
    do c = -grid % n_bnd_cells, -1
      grid % new_c(c) = 0
    end do

    ! Faces step 1/2: inside the domain
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 > 0) then
        if( (grid % comm % cell_proc(c1) .eq. sub) .and.  &
            (grid % comm % cell_proc(c2) .eq. sub) ) then
          nf_sub   = nf_sub + 1
          grid % new_f(s) = nf_sub
        end if
      end if
    end do

    ! Faces step 2/2: on the boundaries
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if( grid % comm % cell_proc(c1) .eq. sub )  then
          nf_sub   = nf_sub + 1
          grid % new_f(s) = nf_sub

          nbc_sub =  nbc_sub + 1       ! increase n. of bnd. cells
          grid % new_c(c2) = -nbc_sub  ! new loc. number of bnd. cell
        end if
      end if 
    end do

    !--------------------!
    !   Create buffers   !
    !--------------------!
    nbf_sub  = 0
    nbfc_sub = 0

    if(verbose) then
      write(fu,'(a30)') '# Number of physical cells:'
      write(fu,'(i8)')  nc_sub
    end if

    do subo = 1, maxval(grid % comm % cell_proc(:))
      if(subo .ne. sub) then
        grid % comm % buff_s_cell(subo) = nbf_sub + 1

        ! Faces inside the domain
        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if(c2  > 0) then
            if( (grid % comm % cell_proc(c1) .eq. sub) .and.  &
                (grid % comm % cell_proc(c2) .eq. subo) ) then
              nbf_sub = nbf_sub + 1                ! increase buffer cell count
              buf_send_ind(nbf_sub) = grid % new_c(c1)  ! buffer send index
              buf_recv_ind(nbf_sub) = c2           ! important for coordinate

              if(grid % new_c(c2) .eq. 0) then
                nbfc_sub = nbfc_sub + 1
                grid % new_c(c2) = nbfc_sub
              end if

              grid % new_f(s) = nf_sub + nbf_sub
            end if
            if( (grid % comm % cell_proc(c2) .eq. sub) .and.  &
                (grid % comm % cell_proc(c1) .eq. subo) ) then
              nbf_sub = nbf_sub + 1                ! increasu buffer cell count
              buf_send_ind(nbf_sub) = grid % new_c(c2)  ! buffer send index
              buf_recv_ind(nbf_sub) = c1           ! important for coordinate

              if(grid % new_c(c1) .eq. 0) then
                nbfc_sub = nbfc_sub + 1
                grid % new_c(c1) = nbfc_sub
              end if

              grid % new_f(s) = nf_sub + nbf_sub
            end if
          end if  ! c2 > 0
        end do    ! through sides

        grid % comm % buff_e_cell(subo) = nbf_sub

        ! Write to buffer file
        if(verbose) then
          write(fu,'(a33)') '#-------------------------------#' 
          write(fu,'(a33)') '#   Conections with subdomain:  #' 
          write(fu,'(a33)') '#-------------------------------#' 
          write(fu,'(i8)')  subo
          write(fu,'(a30)') '# Number of local connections:'
          write(fu,'(i8)')  grid % comm % buff_e_cell(subo) -  &
                            grid % comm % buff_s_cell(subo) + 1
          write(fu,'(a37)') '# Local number in a buffer and index:'
          do b = grid % comm % buff_s_cell(subo),  &
                 grid % comm % buff_e_cell(subo)
            write(fu,'(2i8)') b - grid % comm % buff_s_cell(subo) + 1,  &
                              buf_send_ind(b)
          end do
        end if
      end if

    end do ! for subo

    !-----------!
    !   Nodes   !
    !-----------!
    nn_sub = 0     ! number of cells in subdomain

    ! Initialize node numbers to zero
    do n = 1, grid % n_nodes
      grid % new_n(n) = 0
    end do

    ! Mark nodes for renumbering with -1
    do c = 1, grid % n_cells
      if(grid % comm % cell_proc(c) .eq. sub) then
        do ln = 1, grid % cells_n_nodes(c)
          grid % new_n(grid % cells_n(ln,c)) = -1
        end do
      end if
    end do
    do s = 1, nbf_sub
      do ln = 1, grid % cells_n_nodes(buf_recv_ind(s))
        grid % new_n(grid % cells_n(ln,buf_recv_ind(s))) = -1
      end do
    end do

    ! Renumber marked nodes
    do n = 1, grid % n_nodes
      if(grid % new_n(n) .eq. -1) then
        nn_sub          = nn_sub + 1
        grid % new_n(n) = nn_sub
      end if
    end do

    print '(a,i5,a)', ' #============================================='
    print '(a,i5,a)', ' # Saving subdomain ', sub, ' with:'
    print '(a,i9,a)', ' # ', nc_sub,            ' cells'
    print '(a,i9,a)', ' # ', nn_sub,            ' nodes' 
    print '(a,i9,a)', ' # ', nf_sub + nbf_sub,  ' faces' 
    print '(a,i9,a)', ' # ', nbc_sub,           ' boundary cells' 
    print '(a,i5,a)', ' #---------------------------------------------'

    call Grid_Mod_Save_Cns(grid,         &
                           sub,          &
                           nn_sub,       &
                           nc_sub,       &
                           nf_sub,       &
                           nbc_sub,      &
                           nbf_sub)

    call Grid_Mod_Save_Geo(grid,         &
                           sub,          &
                           nf_sub,       &
                           nbf_sub)

    call Save_Vtu_Cells(grid,       &
                        sub,        &
                        nn_sub,     &
                        nc_sub + nbfc_sub)

    call Save_Vtu_Links(grid,       &
                        sub,        &
                        nn_sub,     &
                        nc_sub,     &
                        nf_sub,     &
                        nbc_sub,    &
                        nbf_sub)

    if(verbose) then
      print '(a)',    ' # Test:'
      print '(a,i9)', ' # nn_sub  :', nn_sub
      print '(a,i9)', ' # nc_sub  :', nc_sub
      print '(a,i9)', ' # nf_sub  :', nf_sub
      print '(a,i9)', ' # nbc_sub :', nbc_sub

      print '(a)',    ' #==============================================='
      print '(a,i9)', ' # Subdomain   ', sub
      print '(a,i9)', ' # Buffer size ', nbf_sub
      do subo = 1, maxval(grid % comm % cell_proc(:))
        if(subo .ne. sub) then
          print '(a,i9,a,3i9)', ' # Connections with ', subo ,' : ',  &
            grid % comm % buff_e_cell(subo) -                         &
            grid % comm % buff_s_cell(subo) + 1,                      &
            nbc_sub + grid % comm % buff_s_cell(subo),                &
            nbc_sub + grid % comm % buff_e_cell(subo)
        end if 
      end do ! for subo
      print '(a)',    ' #-----------------------------------------------'
    end if

  end do   ! through subdomains

  if(verbose) then
    close(fu)
  end if

  !--------------------------------------------------!
  !                                                  !
  !   Save the entire domain with renumbered cells   !
  !                                                  !
  !--------------------------------------------------!
  do n = 1, grid % n_nodes
    grid % new_n(n) = n
  end do
  do c = 1, grid % n_cells
    grid % new_c(c) = 0
  end do
  do s = 1, grid % n_faces
    grid % new_f(s) = 0
  end do

  !------------------------!
  !   Assign new numbers   !
  !------------------------!
  nc_sub = 0     ! number of cells renumbered
  do sub = 1, maxval(grid % comm % cell_proc(:))
    do c = 1, grid % n_cells
      if(grid % comm % cell_proc(c) .eq. sub) then
        nc_sub = nc_sub+1
        grid % new_c(c) = nc_sub
      end if
    end do
  end do

  nf_sub = 0     ! number of sides renumbered
  do sub = 1, maxval(grid % comm % cell_proc(:))
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(grid % comm % cell_proc(c1) .eq. sub) then
        nf_sub = nf_sub+1
        grid % new_f(s) = nf_sub
      end if
    end do
  end do
  print '(a,2i9)', ' # Number of faces: ', grid % n_faces, nf_sub

  !----------------------------------!
  !   Sort what needs to be sorted   !
  !----------------------------------!

  ! Sort cell and face connectivities
  call Grid_Mod_Sort_Cells_By_Index(grid,                &
                                    grid % new_c(1),     &
                                    n=grid % n_cells)
  call Grid_Mod_Sort_Faces_By_Index(grid,                &
                                    grid % new_f(1),     &
                                    n=grid % n_faces)

  ! Sort cell-based values to be plotted
  call Sort_Mod_Int_By_Index(grid % comm % cell_proc(1),  &
                             grid % new_c(1),             &
                             n=grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % wall_dist(1),        &
                              grid % new_c(1),            &
                              n=grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % vol(1),              &
                              grid % new_c(1),            &
                              n=grid % n_cells)

  ! Sort face-cell connectivity
  allocate(side_cell(grid % n_faces, 2))
  do s = 1, grid % n_faces
    side_cell(s,1) = grid % faces_c(1,s)
    side_cell(s,2) = grid % faces_c(2,s)
  end do
  call Sort_Mod_Int_By_Index(side_cell(1,1),      &
                             grid % new_f(1),     &
                             n=grid % n_faces)
  call Sort_Mod_Int_By_Index(side_cell(1,2),      &
                             grid % new_f(1),     &
                             n=grid % n_faces)
  do s = 1, grid % n_faces
    grid % faces_c(1,s) = side_cell(s,1)
    grid % faces_c(2,s) = side_cell(s,2)
  end do
  deallocate(side_cell)

  !----------------------!
  !   And finally save   !
  !----------------------!
  call Save_Vtu_Cells(grid, 0, grid % n_nodes, grid % n_cells)
  call Save_Vtu_Faces(grid)

  end subroutine
