!==============================================================================!
  subroutine Save_Subdomains(Divide, Grid, n_sub)
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!                                                                              !
!   Warning: Unfortunatelly, if you change the order in which cells are        !
!   stored here, it might have an impact on map creation for backup file.      !
!                                                                              !
!   Region-wise, at this point you have:                                       !
!                                                                              !
!   Faces:   1   2   3   4   5   6   7   8   9  10  11  12  13 14  15  16      !
!                                                                              !
!   Cells: -16 -15 -14 -13 -12 -11 -10  -9  -8  -7  -6  -5  -4 -3  -2  -1      !
!                                                                              !
!   Reg:   |<----------- 1 ----------->|<----- 2 ----->|<------ 3 ------>|     !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Divide_Type)  :: Divide
  type(Grid_Type)     :: Grid
  integer, intent(in) :: n_sub  ! number of subdomains
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MARK = -1
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n, s, c1, c2, sub, subo, i_nod
  integer :: ss, sr, nn, reg, bc_indx !exp:, lev
  integer :: nn_sub      ! number of nodes in the subdomain
  integer :: nc_sub      ! number of cells in the subdomain
  integer :: nf_sub      ! number of faces in the subdomain
  integer :: ns_sub      ! number of shadow faces in subdomain
  integer :: nbc_sub     ! number of boundary cells in subdomain
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Divide)
!==============================================================================!

  !-------------------------------!
  !                               !
  !                               !
  !   Browse through subdomains   !
  !                               !
  !                               !
  !-------------------------------!
  do sub = 1, maxval(Grid % Comm % cell_proc(:))

    !------------------!
    !                  !
    !   Cells inside   !
    !                  !
    !------------------!
    nc_sub = 0     ! number of cells in subdomain
    Grid % new_c(1:Grid % n_cells) = 0
    Grid % old_c(1:Grid % n_cells) = 0
    Grid % new_n(1:Grid % n_nodes) = 0

    !--------------------------------------------------------------!
    !   Renumber inside cells the subdomain and mark their nodes   !
    !--------------------------------------------------------------!
    do c = 1, Grid % n_cells
      if(Grid % Comm % cell_proc(c) .eq. sub) then
        nc_sub = nc_sub + 1       ! increase the number of cells in sub.
        Grid % new_c(c) = nc_sub  ! assign new (local) cell number 
        Grid % old_c(nc_sub) = c

        do i_nod = 1, abs(Grid % cells_n_nodes(c))
          Grid % new_n(Grid % cells_n(i_nod,c)) = MARK
        end do
      end if
    end do

    !------------------------------------------!
    !   Also spread node marks to twin nodes   !
    !     (Find them through shadow faces)     !
    !------------------------------------------!
    do subo = 1, maxval(Grid % Comm % cell_proc(:))
      do s = 1, Grid % n_faces
        if(Grid % faces_s(s) > 0) then
          c1 = Grid % faces_c(1, s)
          c2 = Grid % faces_c(2, s)
          if(Grid % Comm % cell_proc(c1) .eq. sub  .and.  &
             Grid % Comm % cell_proc(c2) .eq. subo .or.   &
             Grid % Comm % cell_proc(c2) .eq. sub  .and.  &
             Grid % Comm % cell_proc(c1) .eq. subo) then
            nn = Grid % faces_n_nodes(s)
            ss = Grid % faces_s(s)
            Grid % new_n(Grid % faces_n(1:nn, s )) = MARK
            Grid % new_n(Grid % faces_n(1:nn, ss)) = MARK
          end if
        end if
      end do
    end do

    !-----------------------------!
    !   Inside cells in buffers   !
    !-----------------------------!
    do subo = 1, maxval(Grid % Comm % cell_proc(:))
      if(subo .ne. sub) then

        ! Mark cells in buffer "subo"
        do c = 1, Grid % n_cells
          if(Grid % Comm % cell_proc(c) .eq. subo) then
            n = abs(Grid % cells_n_nodes(c))
            if( any(Grid % new_n(Grid % cells_n(1:n, c)) .eq. MARK) ) then
              Grid % new_c(c) = MARK
            end if
          end if
        end do

        ! Renumber cells in buffer "subo"
        do c = 1, Grid % n_cells
          if(Grid % Comm % cell_proc(c) .eq. subo .and.  &
             Grid % new_c(c) .eq. MARK) then
            nc_sub = nc_sub + 1       ! increase the number of cells in sub.
            Grid % new_c(c) = nc_sub  ! assign new (local) cell number 
            Grid % old_c(nc_sub) = c
          end if
        end do

      end if  ! subo .ne. sub
    end do    ! subo

    !--------------------!
    !                    !
    !   Boundary cells   !
    !                    !
    !--------------------!
    bc_indx = -Grid % n_bnd_cells  ! number of boundary cells in subdomain, ...
                                   ! ... start from the smallest possible
    Grid % new_c(-Grid % n_bnd_cells:-1) = 0
    Grid % old_c(-Grid % n_bnd_cells:-1) = 0

    !---------------------------------------!
    !   Step 1: boundary cells in buffers   !
    !---------------------------------------!
    do subo = 1, maxval(Grid % Comm % cell_proc(:))
      if(subo .ne. sub) then

        do reg = Boundary_Regions()
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            Assert(c2 < 0)
            if( Grid % Comm % cell_proc(c1) .eq. subo .and.  &
                Grid % new_c(c1) .ne. 0)  then
              Grid % new_c(c2) = bc_indx  ! new loc. number of bnd. cell
              bc_indx = bc_indx + 1       ! increase boundary cell counter
            end if
          end do   ! s
        end do     ! region

      end if  ! subo .ne. sub
    end do    ! subo

    !----------------------------------------------!
    !   Step 2: boundary cells in the domain sub   !
    !----------------------------------------------!
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        Assert(c2 < 0)
        if( Grid % Comm % cell_proc(c1) .eq. sub )  then
          Grid % new_c(c2) = bc_indx  ! new loc. number of bnd. cell
          bc_indx = bc_indx + 1       ! increase boundary cell counter
        end if
      end do
    end do

    ! At this point, new boundary cell counters go from -Grid % n_bnd_cells
    ! That is too low, all indices have to be increased in a way that maximum
    ! boundary cell index is -1
    nbc_sub = 0
    do c = -Grid % n_bnd_cells, -1
      if(Grid % new_c(c) .ne. 0)  then
        Grid % new_c(c) = Grid % new_c(c) - bc_indx
        Grid % old_c(Grid % new_c(c)) = c
        nbc_sub = nbc_sub + 1
      end if
    end do

    !-----------!
    !           !
    !   Faces   !
    !           !
    !-----------!
    nf_sub  = 0  ! number of faces in subdomain
    do s = 1, Grid % n_faces + Grid % n_shadows
      Grid % new_f(s) = 0
      Grid % old_f(s) = 0
    end do

    !---------------------------------------------------------!
    !   Step 1: In the domain sub                             !
    !   (It should preserve the order of the original grid.   !
    !    If the faces in the original grid started from       !
    !    boundary faces, it should be the case here too.)     !
    !---------------------------------------------------------!
    do reg = Boundary_And_Inside_Regions()  ! shadows? probably not, at the end
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        ! Both cells are in the domain
        if( Grid % Comm % cell_proc(c1) .eq. sub .and.  &
            Grid % Comm % cell_proc(c2) .eq. sub )  then
          nf_sub = nf_sub + 1
          Grid % new_f(s) = nf_sub
          Grid % old_f(nf_sub) = s
        end if
      end do
    end do

    !-------------------------------------!
    !   Step 2: Inside faces in buffers   !
    !-------------------------------------!
    do subo = 1, maxval(Grid % Comm % cell_proc(:))
      if(subo .ne. sub) then

        ! Faces half in the domain, half in the buffers
        do s = Faces_In_Domain_And_At_Buffers()
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          if( (Grid % Comm % cell_proc(c1) .eq. sub) .and.  &
              (Grid % Comm % cell_proc(c2) .eq. subo) ) then
            nf_sub = nf_sub + 1
            Grid % new_f(s) = nf_sub
            Grid % old_f(nf_sub) = s
          end if
          if( (Grid % Comm % cell_proc(c2) .eq. sub) .and.  &
              (Grid % Comm % cell_proc(c1) .eq. subo) ) then
            nf_sub = nf_sub + 1
            Grid % new_f(s) = nf_sub
            Grid % old_f(nf_sub) = s
          end if

        end do    ! through inside faces (region 0)

      end if  ! subo .ne. sub
    end do    ! subo

    !------------------------------------------------!
    !   Step 3: All the remaining faces in buffers   !
    !------------------------------------------------!
    do reg = Boundary_And_Inside_Regions()  ! shadows? probably not, at the end
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        ! If any of its cells have been marked for saving ...
        if(Grid % new_c(c1) .ne. 0 .or.  &
           Grid % new_c(c2) .ne. 0) then

          ! ... but the face hasn't, do it now!
          if(Grid % new_f(s) .eq. 0) then

            ! Neither of the cells should be in this sub, they are in buffers
            Assert(Grid % Comm % cell_proc(c1) .ne. sub)
            Assert(Grid % Comm % cell_proc(c2) .ne. sub)

            nf_sub = nf_sub + 1
            Grid % new_f(s) = nf_sub
            Grid % old_f(nf_sub) = s
          end if
        end if
      end do
    end do

    !--------------------------!
    !   Step 4: shadow faces   !
    !--------------------------!
    ns_sub  = 0  ! number of shadow faces in subdomain

    ! Browse through shadows only
    do ss = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows

      ! Take real face from the shadow
      sr = Grid % faces_s(ss)

      ! Check if real face was marked for saving
      ! and if it is so, mark also the shadow
      if( Grid % new_f(sr) .ne. 0) then
        ns_sub = ns_sub + 1                 ! increase shadow face counter ...
        Grid % new_f(ss) = nf_sub + ns_sub  ! ... but do not update pointers
        Grid % old_f(nf_sub + ns_sub) = ss  ! from shadow to real and back.
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
    Grid % new_n(1:Grid % n_nodes) = 0

    !---------------------------------------------------!
    !   Mark nodes in cells for renumbering with MARK   !
    !---------------------------------------------------!
    do c = -Grid % n_bnd_cells, Grid % n_cells
      if(Grid % new_c(c) > 0) then
        do i_nod = 1, abs(Grid % cells_n_nodes(c))
          Grid % new_n(Grid % cells_n(i_nod,c)) = MARK
        end do
      end if
    end do

    !---------------------------------------------------!
    !   Mark nodes in faces for renumbering with MARK   !
    !---------------------------------------------------!
    do s = 1, Grid % n_faces + Grid % n_shadows
      if(Grid % new_f(s) > 0) then
        do i_nod = 1, Grid % faces_n_nodes(s)
          Grid % new_n(Grid % faces_n(i_nod,s)) = MARK
        end do
      end if
    end do

    !---------------------------!
    !   Renumber marked nodes   !
    !---------------------------!
    do n = 1, Grid % n_nodes
      if(Grid % new_n(n) .eq. MARK) then
        nn_sub          = nn_sub + 1
        Grid % new_n(n) = nn_sub
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

    call Grid % Save_Cfn((/sub, n_sub/),  &
                           nn_sub,        &
                           nc_sub,        &
                           nf_sub,        &
                           ns_sub,        &   ! number of shadow faces
                           nbc_sub)

    call Grid % Save_Dim((/sub, n_sub/))

    call Grid % Save_Vtu_Cells((/sub, n_sub/),  &
                                 nn_sub,        &
                                 nc_sub)

  end do   ! through subdomains

  end subroutine
