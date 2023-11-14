!==============================================================================!
  subroutine Search_Coordinate_Clusters(Grid, nodal, enforce_uniform,  &
                                        nx, ny, nz)
!------------------------------------------------------------------------------!
!>  This subroutine analyzes a computational grid to identify homogeneous
!>  directions by clustering coordinates. It was developed in collaboration
!>  with ChatGPT and is particularly effective in determining uniformity in
!>  each coordinate direction. Its utility is dependent on two parameters:
!>  MAX_CLUSTERS, which represents the maximum number of clusters in each
!>  direction, and the default round-off factor in Approx_Real.  Enforcement
!>  of uniformity is useful when converting triangular prismatic grids to
!>  polyhedral.  Uniformity of the grid is lost in such cases and this
!>  subroutine can recover it.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization:                                                          !
!     - Prepares arrays to hold cluster data and visibility count.             !
!   * Clustering Process:                                                      !
!     - Iterates over nodes or cells to cluster coordinates, tracking the      !
!       number of occurrences for each distinct value in each direction.       !
!   * Homogeneity Check:                                                       !
!     - Determines if the grid is homogeneous in each direction based on       !
!       cluster uniformity.                                                    !
!   * Optional Grid Modification:                                              !
!     - If enabled, and a homogeneous direction is found, it offers the user   !
!       the option to enforce uniformity in that direction.                    !
!   * Final Report:                                                            !
!     - Provides a summary of the findings, including the number of clusters   !
!       in each direction.                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)               :: Grid   !! grid under consideration
  logical, intent(in)            :: nodal  !! set true to search nodal clusters
  logical, intent(in)            :: enforce_uniform
    !! dictates whether the subroutine should modify the grid
    !! to enforce uniformity in homogeneous directions.
  integer, optional, intent(out) :: nx, ny, nz
    !! output parameter that, if provided, returns the number
    !! of clusters found in the given coordinate direction
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_CLUSTERS = 512
!-----------------------------------[Locals]-----------------------------------!
  integer              :: n, nf, nl, i, clust, c, i_nod
  real                 :: cluster_xyz(MAX_CLUSTERS, 3)
  integer              :: cluster_vis(MAX_CLUSTERS, 3)
  integer              :: cluster_ord(MAX_CLUSTERS)
  real                 :: xyz, length, delta
  integer              :: n_xyz(3) = 0     ! number of clusters in x, y and z
  logical              :: new
  integer, allocatable :: cluster(:)
  character(SL)        :: answer, output
  logical, allocatable :: orphan(:)
!==============================================================================!

  print '(a)', ' #======================================================'
  if(nodal) then
    print '(a)', ' # Searching for nodal coordinate clusters'
  else
    print '(a)', ' # Searching for celular coordinate clusters'
  end if
  print '(a)', ' #------------------------------------------------------'

  ! For some reason, these have to be initialized each time
  cluster_xyz(:,:) = HUGE
  cluster_vis(:,:) = 0
  n_xyz(:) = 0

  if(nodal) then
    nf = 1
    nl = Grid % n_nodes
  else
    nf = 1               ! do not use yet -Grid % n_bnd_cells, if ever
    nl = Grid % n_cells
  end if

  allocate(cluster(nf:nl))
  allocate(orphan (nf:nl))

  ! GMSH might leave some orphan nodes, that ...
  ! ... is nodes which do not belong to any cell
  ! (Is it a problem of Convert itself?)
  if(nodal) then
    orphan(:) = .true.
    do c = 1, Grid % n_cells
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        n = Grid % cells_n(i_nod, c)
        orphan(n) = .false.
      end do
    end do
  ! There should be no orphan cells
  else
    orphan(nf:nl) = .false.
  end if

  !-----------------------------------------!
  !   Browse three cooordinate directions   !
  !-----------------------------------------!
  do i = 1, 3

    ! Browse through all the nodes or cells (including boundary cells)
    do n = nf, nl

      if(.not. orphan(n)) then

        ! Initialize general coordinate xyz
        if(nodal) then
          if(i .eq. 1) xyz = Grid % xn(n)
          if(i .eq. 2) xyz = Grid % yn(n)
          if(i .eq. 3) xyz = Grid % zn(n)
        else
          if(i .eq. 1) xyz = Grid % xc(n)
          if(i .eq. 2) xyz = Grid % yc(n)
          if(i .eq. 3) xyz = Grid % zc(n)
        end if

        ! Find new cluster in each "i" dir
        new = .true.
        do clust = 1, n_xyz(i)
          if( Math % Approx_Real(xyz, cluster_xyz(clust,i)) ) then
            new = .false.  ! nope, you found it
            cluster_vis(clust,i) = cluster_vis(clust,i) + 1
          end if
        end do
        if(new) then
          n_xyz(i) = n_xyz(i) + 1
          if(n_xyz(i) .gt. MAX_CLUSTERS) goto 1   ! this dir goes nowhere
          cluster_xyz(n_xyz(i),i) = xyz
        end if

      end if  ! node not orphan

    end do    ! through nodes in this dir
1   continue
  end do      ! trhough cooridnate directions

  !--------------------------------!
  !   Report what you have found   !
  !--------------------------------!
  write(output, '(a,3i6)') ' # Number of clusters in x, y and z: ', n_xyz(:)
  if(n_xyz(1) > MAX_CLUSTERS) write(output(38:43), '(a6)')  '  high'
  if(n_xyz(2) > MAX_CLUSTERS) write(output(44:49), '(a6)')  '  high'
  if(n_xyz(3) > MAX_CLUSTERS) write(output(50:55), '(a6)')  '  high'
  print '(a)', output

  !-------------------------------------------------------------!
  !   Return the values you have found if parameters are sent   !
  !-------------------------------------------------------------!
  if(present(nx)) then
    nx = n_xyz(1)
    if(nx > MAX_CLUSTERS) nx = -1
  end if

  if(present(ny)) then
    ny = n_xyz(2)
    if(ny > MAX_CLUSTERS) ny = -1
  end if

  if(present(nz)) then
    nz = n_xyz(3)
    if(nz > MAX_CLUSTERS) nz = -1
  end if

  !-----------------------------------------!
  !   Browse three cooordinate directions   !
  !-----------------------------------------!
  do i = 1, 3

    ! If diriection didn't fail
    if(n_xyz(i) .le. MAX_CLUSTERS) then

      ! Hah: found that it's homogeneous
      if( minval(cluster_vis(1:n_xyz(i),i)) .eq.  &
          maxval(cluster_vis(1:n_xyz(i),i)) ) then
        if(i .eq. 1) print '(a,a1)', ' # Grid is homogeneous in x direction'
        if(i .eq. 2) print '(a,a1)', ' # Grid is homogeneous in y direction'
        if(i .eq. 3) print '(a,a1)', ' # Grid is homogeneous in z direction'

        if(nodal .and. enforce_uniform) then
          print '(a)', ' # ... would you like to enforce uniform grid in it?'

          answer = File % Single_Word_From_Keyboard()
          call String % To_Upper_Case(answer)
          if(answer .eq. 'YES') then

            ! Browse through all the nodes to assign them their cluster
            Assert(nf .eq. 1)
            Assert(nl .eq. Grid % n_nodes)
            do n = 1, Grid % n_nodes  ! you get here only in case of nodal
                                      ! meaning nf and nl are superfluous
              if(i .eq. 1) xyz = Grid % xn(n)
              if(i .eq. 2) xyz = Grid % yn(n)
              if(i .eq. 3) xyz = Grid % zn(n)

              do clust = 1, n_xyz(i)
                if( Math % Approx_Real(xyz, cluster_xyz(clust,i)) ) then
                  cluster(n) = clust
                end if
              end do    ! through clusters
            end do      ! through nodes

            ! Order clustered coordinates
            do clust = 1, n_xyz(i)
              cluster_ord(clust) = clust
            end do

            call Sort % Real_Carry_Int(cluster_xyz(1:n_xyz(i),i),  &
                                       cluster_ord(1:n_xyz(i)))

            length = cluster_xyz(n_xyz(i),i) - cluster_xyz(1,i)
            delta  = length / (n_xyz(i) - 1)

            ! Make the coordinates uniform
            do clust = 2, n_xyz(i)-1
              cluster_xyz(clust,i) = cluster_xyz(1,i) + delta * (clust-1)
            end do

            ! Place them in the original order
            call Sort % Real_By_Index(n_xyz(i),                   &
                                      cluster_xyz(1:n_xyz(i),i),  &
                                      cluster_ord(1:n_xyz(i)))

            ! Copy the homogeneized coordinates back to nodes
            do n = 1, Grid % n_nodes
              if(.not. orphan(n)) then
                if(i .eq. 1) Grid % xn(n) = cluster_xyz(cluster(n),i)
                if(i .eq. 2) Grid % yn(n) = cluster_xyz(cluster(n),i)
                if(i .eq. 3) Grid % zn(n) = cluster_xyz(cluster(n),i)
              end if
            end do      ! through nodes

          end if
        end if  ! enforce_uniform
      end if    ! it is homogeneous
    end if      ! direction didn't fail
  end do        ! coordinate directions
  print '(a)', ' #------------------------------------------------------'


  end subroutine
