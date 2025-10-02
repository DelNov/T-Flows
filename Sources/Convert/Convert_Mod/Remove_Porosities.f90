!==============================================================================!
  subroutine Remove_Porosities(Convert, Grid)
!------------------------------------------------------------------------------!
!>  This soubroutine ...
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
!------------------------------[Local parameters]------------------------------!
!-----------------------------------[Locals]-----------------------------------!
  integer              :: dir, c, cnt, bc, i_nod, n
  integer              :: fn(6,4), n_f_nod, f_nod(4), s_nod(4)
  integer              :: c1, c2, s
  logical              :: c1_on, c2_on
  integer, allocatable :: node_on_porosity(:)
  integer, allocatable :: i_work_1(:)
  integer, allocatable :: i_work_2(:,:)
  integer, allocatable :: i_work_3(:,:)
  integer, allocatable :: i_work_4(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Remove_Porosities')

  !--------------------------------!
  !   A couple of initial checks   !
  !--------------------------------!
  Assert(allocated(Grid % old_c))
  Assert(allocated(Grid % new_c))

  !---------------------------------------!
  !                                       !
  !   Phase I: Find cells in porosities   !
  !                                       !
  !---------------------------------------!
  Assert(allocated(Grid % por));

  ! Count the number of cells in porosities (this is for checking only)
  cnt = 0
  do c = 1, Grid % n_cells
    if(Grid % por(c) .gt. 0) cnt = cnt + 1
  end do
  print '(a38,i9)', '# Number of cells in porosities: (2) ', cnt

  !----------------------------!
  !   Mark nodes on porosity   !
  !----------------------------!
  allocate(node_on_porosity(Grid % n_nodes))
  node_on_porosity(:) = 0                    ! a bit of ghost number
  do c = 1, Grid % n_cells
    if(Grid % por(c) .gt. 0) then
      do i_nod = 1, Grid % cells_n_nodes(c)  ! mark node as on building
        n = Grid % cells_n(i_nod, c)
        node_on_porosity(n) = 1              ! ghost number, not great
      end do
    end if
  end do

  !-----------------------------------------------------------!
  !                                                           !
  !   Phase VI: Manage boundary condition names and numbers   !
  !                                                           !
  !-----------------------------------------------------------!

  ! Shift all by one up (to be able to insert building walls as first)
  ! (Check the definition of All_Regions in ../Shared/Browse.h90)
  do bc = Grid % n_bnd_regions + 5, 1, -1
    Grid % region % name(bc+1) = Grid % region % name(bc)
    do c = 1, Grid % n_cells
      do dir = 1, 6
        if(Grid % cells_bnd_region(dir, c) .eq. bc) then
          Grid % cells_bnd_region(dir, c) = bc + 1
        end if
      end do
    end do
  end do

  ! Shift both n_bnd_regions AND n_regions to accomodate
  ! new boundary condition for POROSITY_WALLS.  Although
  ! this works, it is an indication of a flaw in design;
  ! these two variables are dependent of one another,
  ! but can be set independently.
  Grid % n_bnd_regions = Grid % n_bnd_regions + 1
  Grid % n_regions     = Grid % n_regions     + 1

  ! Add first boundary condition for walls
  Grid % region % name(1) = 'POROSITY_WALLS'

  call Grid % Print_Regions_List()

  !------------------------------------------------------!
  !                                                      !
  !   Phase VII: Insert new boundary faces on porosity   !
  !                                                      !
  !------------------------------------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 .gt. 0) then

      c1_on = Grid % por(c1) .eq. 0
      c2_on = Grid % por(c2) .eq. 0

      if(      c1_on .and.       c2_on) cycle
      if(.not. c1_on .and. .not. c2_on) cycle

      if(c1_on .and. .not. c2_on) c = c1
      if(c2_on .and. .not. c1_on) c = c2

      s_nod(1:4) = Grid % faces_n(1:4, s)
      call Sort % Int_Array(s_nod)

      if(Grid % cells_n_nodes(c) .eq. 4) fn = TET
      if(Grid % cells_n_nodes(c) .eq. 5) fn = PYR
      if(Grid % cells_n_nodes(c) .eq. 6) fn = WED
      if(Grid % cells_n_nodes(c) .eq. 8) fn = HEX

      do dir = 1, 6
        if(Grid % cells_bnd_region(dir, c) .eq. 0) then

          n_f_nod    = 0
          f_nod(1:4) = 0

          ! Fill up face nodes
          do i_nod = 1, 4
            if(fn(dir, i_nod) > 0) then
              f_nod(i_nod) = Grid % cells_n(fn(dir, i_nod), c)
              n_f_nod   = n_f_nod + 1
            end if
          end do
          call Sort % Int_Array(f_nod)

          ! Check all nodes from both faves
          if (all(s_nod == f_nod)) then
            Grid % cells_bnd_region(dir, c) = 1
            exit
          end if

        end if  ! cells_bnd_region(dir, c) .eq. 0
      end do    ! dir
    end if      ! c2 .gt .0
  end do        ! through s

  !------------------------------------------------------------!
  !                                                            !
  !   Phase VIII: Store only cells which are not in porosity   !
  !                                                            !
  !------------------------------------------------------------!

  !-----------------------------------------------------!
  !   Store cells' nodes and boundary conditons again   !
  !-----------------------------------------------------!
  allocate(i_work_1(   Grid % n_cells));  i_work_1(:)   = 0
  allocate(i_work_2(8, Grid % n_cells));  i_work_2(:,:) = 0
  allocate(i_work_3(6, Grid % n_cells));  i_work_3(:,:) = 0
  allocate(i_work_4(   Grid % n_cells));  i_work_4(:)   = 0

  do c = 1, Grid % n_cells
    i_work_1(c)      = Grid % cells_n_nodes(c)
    i_work_2(1:8, c) = Grid % cells_n(1:8, c)
    i_work_3(1:6, c) = Grid % cells_bnd_region(1:6, c)
    i_work_4(c)      = Grid % por(c)
  end do

  !------------------------------------------------------!
  !   Renumber the cells without the cells in porosity   !
  !------------------------------------------------------!
  Grid % new_c(:) = 0
  cnt             = 0
  do c = 1, Grid % n_cells
    if(Grid % por(c) == 0) then
      cnt = cnt + 1
      Grid % new_c(c) = cnt
    end if
  end do
  print '(a38,i9)', '# Old number of cells:               ', Grid % n_cells
  print '(a38,i9)', '# Number of cells without buildings: ', cnt
  print '(a38,i9)', '# Number of excluded cells:          ', Grid % n_cells-cnt

  !--------------------------------------!
  !   Information on cell connectivity   !
  !--------------------------------------!
  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. 0) then
      Grid % cells_n_nodes        (Grid % new_c(c)) = i_work_1(c)
      Grid % cells_n         (1:8, Grid % new_c(c)) = i_work_2(1:8, c)
      Grid % cells_bnd_region(1:6, Grid % new_c(c)) = i_work_3(1:6, c)
      Grid % por                  (Grid % new_c(c)) = i_work_4(c)
    end if
  end do
  Grid % n_cells = cnt

  call Profiler % Stop('Remove_Porosities')

  end subroutine
