!==============================================================================!
  subroutine Save_Subdomains(grid, n_buff_layers)
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!                                                                              !
!   Warning: Unfortunatelly, if you change the order in which cells are        !
!   stored here, it might have an impact on map creation for backup file.      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
  use Grid_Mod,      only: Grid_Type,                     &
                           Grid_Mod_Sort_Cells_By_Index,  &
                           Grid_Mod_Sort_Faces_By_Index,  &
                           Grid_Mod_Save_Cfn,             &
                           Grid_Mod_Save_Dim
  use Sort_Mod       ! it's a collection of subroutines, no need for "only"
  use Save_Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n_buff_layers  ! number of buffer layers
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n, s, c1, c2, sub, subo, i_nod, i_fac, ss, sr, nn !exp:, lev
  integer :: nn_sub      ! number of nodes in the subdomain
  integer :: nc_sub      ! number of cells in the subdomain
  integer :: nf_sub      ! number of faces in the subdomain
  integer :: ns_sub      ! number of shadow faces in subdomain
  integer :: nbc_sub     ! number of boundary cells in subdomain
!==============================================================================!

  !-------------------------------!
  !                               !
  !                               !
  !   Browse through subdomains   !
  !                               !
  !                               !
  !-------------------------------!
  do sub = 1, maxval(grid % comm % cell_proc(:))

    !------------------!
    !                  !
    !   Cells inside   !
    !                  !
    !------------------!
    nc_sub = 0     ! number of cells in subdomain
    grid % new_c(1:grid % n_cells) = 0
    grid % old_c(1:grid % n_cells) = 0
    grid % new_n(1:grid % n_nodes) = 0

    !--------------------------------------------------------------!
    !   Renumber inside cells the subdomain and mark their nodes   !
    !--------------------------------------------------------------!
    do c = 1, grid % n_cells
      if(grid % comm % cell_proc(c) .eq. sub) then
        nc_sub = nc_sub + 1       ! increase the number of cells in sub.
        grid % new_c(c) = nc_sub  ! assign new (local) cell number 
        grid % old_c(nc_sub) = c

        do i_nod = 1, abs(grid % cells_n_nodes(c))
          grid % new_n(grid % cells_n(i_nod,c)) = -1
        end do
      end if
    end do

    !------------------------------------------!
    !   Also spread node marks to twin nodes   !
    !------------------------------------------!
    do subo = 1, maxval(grid % comm % cell_proc(:))
      do s = 1, grid % n_faces
        if(grid % faces_s(s) > 0) then
          c1 = grid % faces_c(1, s)
          c2 = grid % faces_c(2, s)
          if(grid % comm % cell_proc(c1) .eq. sub  .and.  &
             grid % comm % cell_proc(c2) .eq. subo .or.   &
             grid % comm % cell_proc(c2) .eq. sub  .and.  &
             grid % comm % cell_proc(c1) .eq. subo) then
            nn = grid % faces_n_nodes(s)
            ss = grid % faces_s(s)
            grid % new_n(grid % faces_n(1:nn, s )) = -1
            grid % new_n(grid % faces_n(1:nn, ss)) = -1
          end if
        end if
      end do
    end do

    !-----------------------------!
    !   Inside cells in buffers   !
    !-----------------------------!
    do subo = 1, maxval(grid % comm % cell_proc(:))
      if(subo .ne. sub) then

        ! Mark cells in buffer "subo"
        do c = 1, grid % n_cells
          if(grid % comm % cell_proc(c) .eq. subo) then
            n = abs(grid % cells_n_nodes(c))
            if( any(grid % new_n(grid % cells_n(1:n, c)) .eq. -1) ) then
              grid % new_c(c) = -1
            end if
          end if
        end do

!exp:        ! Browse through deeper levels of buffers
!exp:        do lev = 2, n_buff_layers
!exp:
!exp:          ! Mark nodes on this level ...
!exp:          do c = 1, grid % n_cells
!exp:            if(grid % new_c(c) .eq. -1) then
!exp:              do i_nod = 1, abs(grid % cells_n_nodes(c))
!exp:                grid % new_n(grid % cells_n(i_nod,c)) = -1
!exp:              end do
!exp:            end if
!exp:          end do
!exp:
!exp:          ! ... and then also the cells
!exp:          do c = 1, grid % n_cells
!exp:            if(grid % comm % cell_proc(c) .eq. subo) then
!exp:              n = abs(grid % cells_n_nodes(c))
!exp:              if( any(grid % new_n(grid % cells_n(1:n, c)) .eq. -1) ) then
!exp:                grid % new_c(c) = -1
!exp:              end if
!exp:            end if
!exp:          end do
!exp:        end do

        ! Renumber cells in buffer "subo"
        do c = 1, grid % n_cells
          if(grid % comm % cell_proc(c) .eq. subo .and.  &
             grid % new_c(c) .eq. -1) then
            nc_sub = nc_sub + 1       ! increase the number of cells in sub.
            grid % new_c(c) = nc_sub  ! assign new (local) cell number 
            grid % old_c(nc_sub) = c
          end if
        end do

      end if  ! subo .ne. sub
    end do    ! subo

    !--------------------!
    !                    !
    !   Boundary cells   !
    !                    !
    !--------------------!
    nbc_sub = 0     ! number of boundary cells in subdomain
    grid % new_c(-grid % n_bnd_cells:-1) = 0
    grid % old_c(-grid % n_bnd_cells:-1) = 0

    !---------------------------------------!
    !   Step 1: boundary cells in buffers   !
    !---------------------------------------!
    do subo = 1, maxval(grid % comm % cell_proc(:))
      if(subo .ne. sub) then

        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if(c2 < 0) then
            if( grid % comm % cell_proc(c1) .eq. subo .and.  &
                grid % new_c(c1) .ne. 0)  then
              nbc_sub = nbc_sub + 1        ! increase n. of bnd. cells
              grid % new_c(c2) = -nbc_sub  ! new loc. number of bnd. cell
              grid % old_c(-nbc_sub) = c2
            end if
          end if
        end do

      end if  ! subo .ne. sub
    end do    ! subo

    !----------------------------------------------!
    !   Step 2: boundary cells in the domain sub   !
    !----------------------------------------------!
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if( grid % comm % cell_proc(c1) .eq. sub )  then
          nbc_sub = nbc_sub + 1        ! increase n. of bnd. cells
          grid % new_c(c2) = -nbc_sub  ! new loc. number of bnd. cell
          grid % old_c(-nbc_sub) = c2
        end if
      end if
    end do

    !-----------!
    !           !
    !   Faces   !
    !           !
    !-----------!
    nf_sub  = 0  ! number of faces in subdomain
    do s = 1, grid % n_faces + grid % n_shadows
      grid % new_f(s) = 0
      grid % old_f(s) = 0
    end do

    !---------------------------------------------------------!
    !   Step 1: In the domain sub                             !
    !   (It should preserve the order of the original grid.   !
    !    If the faces in the original grid started from       !
    !    boundary faces, it should be the case here too.)     !
    !---------------------------------------------------------!
    do s = 1, grid % n_faces ! + n_shadows???
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Both cells are in the domain
      if( grid % comm % cell_proc(c1) .eq. sub .and.  &
          grid % comm % cell_proc(c2) .eq. sub )  then
        nf_sub = nf_sub + 1
        grid % new_f(s) = nf_sub
        grid % old_f(nf_sub) = s
      end if
    end do

    !-------------------------------------!
    !   Step 2: Inside faces in buffers   !
    !-------------------------------------!
    do subo = 1, maxval(grid % comm % cell_proc(:))
      if(subo .ne. sub) then

        ! Faces half in the domain, half in the buffers
        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if(c2 > 0) then
            if( (grid % comm % cell_proc(c1) .eq. sub) .and.  &
                (grid % comm % cell_proc(c2) .eq. subo) ) then
              nf_sub = nf_sub + 1
              grid % new_f(s) = nf_sub
              grid % old_f(nf_sub) = s
            end if
            if( (grid % comm % cell_proc(c2) .eq. sub) .and.  &
                (grid % comm % cell_proc(c1) .eq. subo) ) then
              nf_sub = nf_sub + 1
              grid % new_f(s) = nf_sub
              grid % old_f(nf_sub) = s
            end if
          end if  ! c2 > 0
        end do    ! through faces

      end if  ! subo .ne. sub
    end do    ! subo

    !------------------------------------------------!
    !   Step 3: All the remaining faces in buffers   !
    !------------------------------------------------!
    do c = -grid % n_bnd_cells, grid % n_cells

      ! If cell has been marked to be saved with this sub, no matter
      ! if it is inside or at the buffer (although it would be strange
      ! if inside faces haven't been marked already above)
      if(grid % new_c(c) .ne. 0) then

        do i_fac = 1, grid % cells_n_faces(c)
          s = grid % cells_f(i_fac, c)

          ! If the face hasn't been marked allready, do it now!
          if(grid % new_f(s) .eq. 0) then
            nf_sub = nf_sub + 1
            grid % new_f(s) = nf_sub
            grid % old_f(nf_sub) = s
          end if
        end do  ! i_fac

      end if  !  new_c(c) .ne. 0
    end do    !  c

    !--------------------------!
    !   Step 4: shadow faces   !
    !--------------------------!
    ns_sub  = 0  ! number of shadow faces in subdomain

    ! Browse through shadows only
    do ss = grid % n_faces + 1, grid % n_faces + grid % n_shadows

      ! Take real face from the shadow
      sr = grid % faces_s(ss)

      ! Check if real face was marked for saving
      ! and if it is so, mark also the shadow
      if( grid % new_f(sr) .ne. 0) then
        ns_sub = ns_sub + 1                 ! increase shadow face counter ...
        grid % new_f(ss) = nf_sub + ns_sub  ! ... but do not update pointers
        grid % old_f(nf_sub + ns_sub) = ss  ! from shadow to real and back.
      end if

    end do    ! through faces

    !-----------!
    !           !
    !   Nodes   !
    !           !
    !-----------!
    nn_sub = 0     ! number of cells in subdomain

    !-------------------------------------!
    !   Initialize node numbers to zero   !
    !-------------------------------------!
    grid % new_n(1:grid % n_nodes) = 0

    !-------------------------------------------------!
    !   Mark nodes in cells for renumbering with -1   !
    !-------------------------------------------------!
    do c = -grid % n_bnd_cells, grid % n_cells
      if(grid % new_c(c) > 0) then
        do i_nod = 1, abs(grid % cells_n_nodes(c))
          grid % new_n(grid % cells_n(i_nod,c)) = -1
        end do
      end if
    end do

    !-------------------------------------------------!
    !   Mark nodes in faces for renumbering with -1   !
    !-------------------------------------------------!
    do s = grid % n_faces + 1, grid % n_faces + grid % n_shadows
      if(grid % new_f(s) > 0) then
        do i_nod = 1, grid % faces_n_nodes(s)
          grid % new_n(grid % faces_n(i_nod,s)) = -1
        end do
      end if
    end do

    !---------------------------!
    !   Renumber marked nodes   !
    !---------------------------!
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
    print '(a,i9,a)', ' # ', nf_sub,            ' faces'
    print '(a,i9,a)', ' # ', ns_sub,            ' shadow faces'
    print '(a,i9,a)', ' # ', nbc_sub,           ' boundary cells'
    print '(a,i5,a)', ' #---------------------------------------------'

    call Grid_Mod_Save_Cfn(grid,         &
                           sub,          &
                           nn_sub,       &
                           nc_sub,       &
                           nf_sub,       &
                           ns_sub,       &   ! number of shadow faces
                           nbc_sub)

    call Grid_Mod_Save_Dim(grid,         &
                           sub)

    call Save_Vtu_Cells(grid,       &
                        sub,        &
                        nn_sub,     &
                        nc_sub)

  end do   ! through subdomains

  end subroutine
