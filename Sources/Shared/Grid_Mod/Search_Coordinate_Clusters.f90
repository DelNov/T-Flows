!==============================================================================!
  subroutine Search_Coordinate_Clusters(Grid, enforce_uniform)
!------------------------------------------------------------------------------!
!   This subroutine searches for homogeneous directions in a grid.  It is a    !
!   result of discussion with ChatGPT, which did come with a good idea to      !
!   cluster coorindates and in that way check uniformity of a grid in each     !
!   particular direction.                                                      !
!                                                                              !
!   Although elegant and sleak, the subroutine does depend on two parameters;  !
!   namely the MAX_CLUSTERS, maximum number of clusers in each coordinate      !
!   direction (which would be analogous to number of cells in each direction   !
!   if grid was Cartesian) and the default round-off factor in Approx_Real,    !
!   which is set to 1.0e-9 by default.  Although these factors work in 2023,   !
!   they migh have to be re-visited in the future if CPU power goes up a lot   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  logical, intent(in) :: enforce_uniform
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_CLUSTERS = 512
!-----------------------------------[Locals]-----------------------------------!
  integer              :: n, i, clust, C, I_NOD
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
  print '(a)', ' # Searching for nodal coordinate clusters'
  print '(a)', ' #------------------------------------------------------'

  ! For some reason, these have to be initialized each time
  cluster_xyz(:,:) = HUGE
  cluster_vis(:,:) = 0
  n_xyz(:) = 0

  allocate(cluster(Grid % n_nodes))
  allocate(orphan (Grid % n_nodes))
  orphan(:) = .true.

  ! GMSH might leave some orphan nodes, that ...
  ! ... is nodes which do not belong to any cell
  ! (Is it a problem of Convert itself?)
  do c = 1, Grid % n_cells
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      orphan(n) = .false.
    end do
  end do

  !-----------------------------------------!
  !   Browse three cooordinate directions   !
  !-----------------------------------------!
  do i = 1, 3

    ! Browse through all the nodes
    do n = 1, Grid % n_nodes

      if(.not. orphan(n)) then

        ! Initialize general coordinate xyz
        if(i .eq. 1) xyz = Grid % xn(n)
        if(i .eq. 2) xyz = Grid % yn(n)
        if(i .eq. 3) xyz = Grid % zn(n)

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

        if(enforce_uniform) then
          print '(a)', ' # ... would you like to enforce uniform grid in it?'

          answer = File % Single_Word_From_Keyboard()
          call String % To_Upper_Case(answer)
          if(answer .eq. 'YES') then

            ! Browse through all the nodes to assign them their cluster
            do n = 1, Grid % n_nodes

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
