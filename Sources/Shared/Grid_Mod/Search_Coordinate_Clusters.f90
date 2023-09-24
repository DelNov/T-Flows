!==============================================================================!
  subroutine Search_Coordinate_Clusters(Grid)
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
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer      :: n, i, cluster
  real         :: cluster_xyz(3, MAX_CLUSTERS) = HUGE
  integer      :: cluster_vis(3, MAX_CLUSTERS) = 0
  real         :: xyz                         ! x, y, or z coordinate
  integer      :: n_xyz(3) = 0                ! number of clusters in x, y and z
  logical      :: new
!==============================================================================!

  print '(a)', ' #======================================================'
  print '(a)', ' # Searching for nodal coordinate clusters'
  print '(a)', ' #------------------------------------------------------'

  !-----------------------------------------!
  !   Browse three cooordinate directions   !
  !-----------------------------------------!
  do i = 1, 3

    ! Browse through all the nodes
    do n = 1, Grid % n_nodes

      ! Initialize general coordinate xyz
      if(i .eq. 1) xyz = Grid % xn(n)
      if(i .eq. 2) xyz = Grid % yn(n)
      if(i .eq. 3) xyz = Grid % zn(n)

      ! Find new cluster in each "i" dir
      new = .true.
      do cluster = 1, n_xyz(i)
        if( Math % Approx_Real(xyz, cluster_xyz(i, cluster)) ) then
          new = .false.  ! nope, you found it
          cluster_vis(i, cluster) = cluster_vis(i, cluster) + 1
        end if
      end do
      if(new) then
        n_xyz(i) = n_xyz(i) + 1
        if(n_xyz(i) .gt. MAX_CLUSTERS) goto 1   ! this dir goes nowhere
        cluster_xyz(i,n_xyz(i)) = xyz
      end if

    end do  ! through nodes in this dir
1   continue
  end do

  !--------------------------------!
  !   Report what you have found   !
  !--------------------------------!
  print '(a,3i6)',   ' # Number of clusters in x, y and z: ', n_xyz(:)

  do i = 1, 3
    if(n_xyz(i) .le. MAX_CLUSTERS) then
      if( minval(cluster_vis(i, 1:n_xyz(i))) .eq.  &
          maxval(cluster_vis(i, 1:n_xyz(i))) ) then
        if(i .eq. 1) print '(a,a1)', ' # Grid seems to be homogeneous in x dir'
        if(i .eq. 2) print '(a,a1)', ' # Grid seems to be homogeneous in y dir'
        if(i .eq. 3) print '(a,a1)', ' # Grid seems to be homogeneous in z dir'
      end if
    end if
  end do

  end subroutine
