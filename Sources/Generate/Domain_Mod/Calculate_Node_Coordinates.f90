!==============================================================================!
  subroutine Calculate_Node_Coordinates(Dom, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine, used only from Generate program, is responsible for
!>  calculating the coordinates of nodes within each block of the computational
!>  domain and transferring this data to a grid object.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initialization: It first initializes the node coordinates to a specific  !
!     value (PETA) to serve as an indicator of whether a node's coordinates    !
!     have been calculated.                                                    !
!   * Block-by-block processing: For each block in the domain, it calculates   !
!     the coordinates of the nodes based on the block's resolution and corner  !
!     points.                                                                  !
!   * Processing lines within blocks: The subroutine then iterates over each   !
!     line within the domain, checking if the line is part of the current      !
!     block and then calculating node coordinates along the line based on      !
!     either direct specifications (point by point) or a weighting factor.     !
!   * Distributing nodes on faces and volumes: It distributes nodes along the  !
!     faces and within the volumes of each block, considering different        !
!     weighting factors for node clustering.                                   !
!   * Control volume nodes and neighbors: Finally, the subroutine calculates   !
!     and assigns nodes and their neighbors to the control volumes within      !
!     each block, updating the overall number of nodes and cells in the grid.  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type) :: Dom   !! domain in which the grid is being generated
  type(Grid_Type)    :: Grid  !! grid being generated
!-----------------------------------[Locals]-----------------------------------!
  integer :: fc, b, bl, i, j, k, n, c, ig
  integer :: l, l1, l2
  integer :: is, js, ks, ie, je, ke
  integer :: ni, nj, nk, ci, cj, ck
  integer :: trans(3,2)
!==============================================================================!

  ! Setting node_positioned to false here serves later as indicator if
  ! coordinates this particular node has been calculated or not
  ! (For example in Domain_Mod_Distribute_Nodes Domain_Mod_Laplace)
  do n = 1, Grid % max_n_nodes
    Grid % node_positioned(n) = .false.
  end do

  !------------------------------------!
  !   Calculate the node coordinates   !
  !------------------------------------!
  Grid % n_nodes = 0  ! initialize n.o.n.
  Grid % n_cells = 0  ! initialize n.o.v.

  do b = 1, size(Dom % blocks)

    print '(a38,i7)', '# Generating block:                  ', b
    ni = Dom % blocks(b) % resolutions(1)
    nj = Dom % blocks(b) % resolutions(2)
    nk = Dom % blocks(b) % resolutions(3)

    ! ( 1 )
    n = Grid % n_nodes + ( 1-1)*ni*nj + ( 1-1)*ni +  1
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(1)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(1)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(1)) % z

    ! ( 2 )
    n = Grid % n_nodes + ( 1-1)*ni*nj + ( 1-1)*ni + ni
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(2)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(2)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(2)) % z

    ! ( 3 )
    n = Grid % n_nodes + ( 1-1)*ni*nj + (nj-1)*ni +  1
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(3)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(3)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(3)) % z

    ! ( 4 )
    n = Grid % n_nodes + ( 1-1)*ni*nj + (nj-1)*ni + ni
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(4)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(4)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(4)) % z

    ! ( 5 ) !
    n = Grid % n_nodes + (nk-1)*ni*nj + ( 1-1)*ni +  1
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(5)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(5)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(5)) % z

    ! ( 6 ) !
    n = Grid % n_nodes + (nk-1)*ni*nj + ( 1-1)*ni + ni
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(6)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(6)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(6)) % z

    ! ( 7 ) !
    n = Grid % n_nodes + (nk-1)*ni*nj + (nj-1)*ni +  1
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(7)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(7)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(7)) % z

    ! ( 8 ) !
    n = Grid % n_nodes + (nk-1)*ni*nj + (nj-1)*ni + ni
    Grid % xn(n) = Dom % points(Dom % blocks(b) % corners(8)) % x
    Grid % yn(n) = Dom % points(Dom % blocks(b) % corners(8)) % y
    Grid % zn(n) = Dom % points(Dom % blocks(b) % corners(8)) % z

    !------------------------------!
    !   First on the Dom % lines   !
    !    defined point by point    !
    !------------------------------!
    do l=1, size(Dom % lines)

      bl = Dom % Is_Line_in_Block(Dom % lines(l) % points(1),  &
                                  Dom % lines(l) % points(2),  &
                                  b)

      if(bl .eq. b) then

        do n = 1, 8
          if(      Dom % lines(l)  % points(1)  &
              .eq. Dom % blocks(b) % corners(n) ) l1=n
          if(      Dom % lines(l)  % points(2)  &
              .eq. Dom % blocks(b) % corners(n) ) l2=n
        end do

        ! Line is defined in the +i direction
        if     ( (l2-l1) .eq. +1 ) then
          trans(1,1) = 0
          trans(1,2) =+1
          trans(2,2) = 0
          trans(3,2) = 0
          if( (l1 .eq. 1).or.(l1 .eq. 5) ) then
            trans(2,1) =1
          else
            trans(2,1) =Dom % blocks(b) % resolutions(2)
          end if
          if( (l1 .eq. 1).or.(l1 .eq. 3) ) then
            trans(3,1) =1
          else
            trans(3,1) =Dom % blocks(b) % resolutions(3)
          end if

        ! Line is defined in the -i direction
        else if( (l2-l1) .eq. -1 ) then
          trans(1,1) =Dom % blocks(b) % resolutions(1)+1   ! ni from block + 1
          trans(1,2) =-1
          trans(2,2) = 0
          trans(3,2) = 0
          if( (l1 .eq. 2).or.(l1 .eq. 6) ) then
            trans(2,1) =1
          else
            trans(2,1) =Dom % blocks(b) % resolutions(2)
          end if
          if( (l1 .eq. 2).or.(l1 .eq. 4) ) then
            trans(3,1) =1
          else
            trans(3,1) =Dom % blocks(b) % resolutions(3)
          end if

        ! Line is defined in the +j direction
        else if( (l2-l1) .eq. +2 ) then
          trans(2,1) = 0
          trans(2,2) =+1
          trans(1,2) = 0
          trans(3,2) = 0
          if( (l1 .eq. 1).or.(l1 .eq. 5) ) then
            trans(1,1) =1
          else
            trans(1,1) =Dom % blocks(b) % resolutions(1)
          end if
          if( (l1 .eq. 1).or.(l1 .eq. 2) ) then
            trans(3,1) =1
          else
            trans(3,1) =Dom % blocks(b) % resolutions(3)
          end if

        ! Line is defined in the -j direction
        else if( (l2-l1) .eq. -2 ) then
          trans(2,1) =Dom % blocks(b) % resolutions(2)+1   ! nj from block + 1
          trans(2,2) =-1
          trans(1,2) = 0
          trans(3,2) = 0
          if( (l1 .eq. 3).or.(l1 .eq. 7) ) then
            trans(1,1) =1
          else
            trans(1,1) =Dom % blocks(b) % resolutions(1)
          end if
          if( (l1 .eq. 3).or.(l1 .eq. 4) ) then
            trans(3,1) =1
          else
            trans(3,1) =Dom % blocks(b) % resolutions(3)
          end if

        ! Line is defined in the +k direction
        else if( (l2-l1) .eq. +4 ) then
          trans(3,1) = 0
          trans(3,2) =+1
          trans(1,2) = 0
          trans(2,2) = 0
          if( (l1 .eq. 1).or.(l1 .eq. 3) ) then
            trans(1,1) =1
          else
            trans(1,1) =Dom % blocks(b) % resolutions(1)
          end if
          if( (l1 .eq. 1).or.(l1 .eq. 2) ) then
            trans(2,1) =1
          else
            trans(2,1) =Dom % blocks(b) % resolutions(2)
          end if

        ! Line is defined in the -k direction
        else if( (l2-l1) .eq. -4 ) then
          trans(3,1) =Dom % blocks(b) % resolutions(3) + 1  ! nk from block + 1
          trans(3,2) =-1
          trans(1,2) = 0
          trans(2,2) = 0
          if( (l1 .eq. 5).or.(l1 .eq. 7) ) then
            trans(1,1) =1
          else
            trans(1,1) =Dom % blocks(b) % resolutions(1)
          end if
          if( (l1 .eq. 5).or.(l1 .eq. 6) ) then
            trans(2,1) =1
          else
            trans(2,1) =Dom % blocks(b) % resolutions(2)
          end if

        end if ! l1-l2

        ! Line is defined point by point
        if( Math % Approx_Real( Dom % lines(l) % weight, 0.0) ) then
          print *, '# Line: ', l
          print *, '# l1= ', l1
          print *, '# l2= ', l2
          do ig=1,Dom % lines(l) % resolution
            i=trans(1,1)+trans(1,2)*ig
            j=trans(2,1)+trans(2,2)*ig
            k=trans(3,1)+trans(3,2)*ig

            n = Grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i
            Grid % xn(n) = Dom % lines(l) % x(ig)
            Grid % yn(n) = Dom % lines(l) % y(ig)
            Grid % zn(n) = Dom % lines(l) % z(ig)
          end do

        ! Line is defined with a weight factor
        else
          is=trans(1,1)+trans(1,2)
          js=trans(2,1)+trans(2,2)
          ks=trans(3,1)+trans(3,2)
          ie=trans(1,1)+trans(1,2)*Dom % lines(l) % resolution
          je=trans(2,1)+trans(2,2)*Dom % lines(l) % resolution
          ke=trans(3,1)+trans(3,2)*Dom % lines(l) % resolution
          call Dom % Distribute_Nodes(Grid,                  &
                                b, Dom % lines(l) % weight,  &
                                is, js, ks, ie, je, ke)
        end if

      end if ! if the block contains

    end do ! for the Dom % lines

    !-----------!
    !   Lines   !
    !-----------!
    do k=1,nk,nk-1
      do j=1,nj,nj-1
        call Dom % Distribute_Nodes(Grid,  &
                              b, Dom % blocks(b) % weights(1), 1,j,k,ni,j,k)
      end do
    end do

    do k=1,nk,nk-1
      do i=1,ni,ni-1
        call Dom % Distribute_Nodes(Grid,  &
                              b, Dom % blocks(b) % weights(2), i,1,k,i,nj,k)
      end do
    end do

    do j=1,nj,nj-1
      do i=1,ni,ni-1
        call Dom % Distribute_Nodes(Grid,  &
                              b, Dom % blocks(b) % weights(3), i,j,1,i,j,nk)
      end do
    end do

    !---------------------------------------------------------------!
    !   Surfaces...                                                 !
    !                                                               !
    !   I think this is the propper way to calculate surfaces:      !
    !   it spans the Dom % lines in the direction of higher weigh   !
    !---------------------------------------------------------------!

    ! I (k=1)
    fc = 1   ! face index
    k = 1
    if( .not. Math % Approx_Real(  &
              Dom % blocks(b) % face_weights(fc,1),1.0 ) ) then
      do j=1,nj
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    else ! Dom % lines in the j direction
      do i=1,ni
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    end if

    ! VI (k=nk)
    fc = 6   ! face index
    k = nk
    if( .not. Math % Approx_Real(  &
              Dom % blocks(b) % face_weights(fc,1),1.0 ) ) then
     do j=1,nj
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    else ! Dom % lines in the j direction
      do i=1,ni
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    end if

    ! V (i=1)
    fc = 5   ! face index
    i = 1
    if( .not. Math % Approx_Real(  &
              Dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do j=1,nj
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,3),  &
                              i,j,1,i,j,nk)
      end do
    else ! Dom % lines in the j direction
      do k=1,nk
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    end if

    ! III (i=ni)
    fc = 3   ! face index
    i = ni
    if( .not. Math % Approx_Real(  &
              Dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do j=1,nj
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,3),  &
                              i,j,1,i,j,nk)
      end do
    else ! Dom % lines in the j direction
      do k=1,nk
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    end if

    ! II (j=1)
    fc = 2   ! face index
    j = 1
    if( .not. Math % Approx_Real(  &
              Dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do i=1,ni
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,3),  &
                              i,j,1,i,j,nk)
      end do
    else ! Dom % lines in the i direction
      do k=1,nk
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    end if

    ! IV (j=nj)
    fc = 4   ! face index
    j = nj
    if( .not. Math % Approx_Real(  &
              Dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do i=1,ni
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,3),  &
                              i,j,1,i,j,nk)
      end do
    else ! Dom % lines in the i direction
      do k=1,nk
        call Dom % Distribute_Nodes(Grid,                               &
                              b, Dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    end if

    !-------------!
    !   Volumes   !
    !-------------!
    if( .not. Math % Approx_Real(  &
              Dom % blocks(b) % weights(3), 1.0 ) ) then
      do i=1,ni
        do j=1,nj
          call Dom % Distribute_Nodes(Grid,                           &
                                b, Dom % blocks(b) % weights(3), i,j,1,i,j,nk)
        end do
      end do
    else if( .not. Math % Approx_Real(  &
                   Dom % blocks(b) % weights(1), 1.0 ) ) then
      do k=1,nk
        do j=1,nj
          call Dom % Distribute_Nodes(Grid,                           &
                                b, Dom % blocks(b) % weights(1), 1,j,k,ni,j,k)
        end do
      end do
    else if( .not. Math % Approx_Real(  &
                   Dom % blocks(b) % weights(2), 1.0 ) ) then
      do k=1,nk
        do i=1,ni
          call Dom % Distribute_Nodes(Grid,                           &
                                b, Dom % blocks(b) % weights(2), i,1,k,i,nj,k)
        end do
      end do
    else

      do i=1,ni
        do j=1,nj
          do k=1,nk
            n = Grid % n_nodes+(k-1)*ni*nj + (j-1)*ni + i
            call Dom % Laplace(Grid, b, i, j, k,                 &
                               ONE_THIRD, ONE_THIRD, ONE_THIRD,  &
                               ONE_THIRD, ONE_THIRD, ONE_THIRD,  &
                               ONE_THIRD, ONE_THIRD, ONE_THIRD)
          end do
        end do
      end do
    end if

    !-----------------------------------------!
    !   Set the control volume nodes (CellN)  !
    !    and  the control volume neighbours   !
    !-----------------------------------------!
    ci = ni-1
    cj = nj-1
    ck = nk-1

    do k=1,ck
      do j=1,cj
        do i=1,ci
          c = Grid % n_cells + (k-1)*ci*cj + (j-1)*ci + i ! cell
          n = Grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i ! 1st node

          ! Nodes
          call Enlarge % Matrix_Int(Grid % cells_n, i=(/1,8/))
          Grid % cells_n(1,c) = n
          Grid % cells_n(2,c) = n+1
          Grid % cells_n(3,c) = n+ni
          Grid % cells_n(4,c) = n+ni+1
          Grid % cells_n(5,c) = Grid % cells_n(1,c)+ni*nj
          Grid % cells_n(6,c) = Grid % cells_n(2,c)+ni*nj
          Grid % cells_n(7,c) = Grid % cells_n(3,c)+ni*nj
          Grid % cells_n(8,c) = Grid % cells_n(4,c)+ni*nj

          ! Neighbours
          call Enlarge % Matrix_Int(Grid % cells_c, i=(/1,6/))
          Grid % cells_c(1,c) = c-ci*cj
          Grid % cells_c(2,c) = c-ci
          Grid % cells_c(3,c) = c+1
          Grid % cells_c(4,c) = c+ci
          Grid % cells_c(5,c) = c-1
          Grid % cells_c(6,c) = c+ci*cj

          ! This value (-1) is also the default boundary marker
          if(i .eq. 1)  Grid % cells_c(5,c) =-1
          if(i .eq. ci) Grid % cells_c(3,c) =-1
          if(j .eq. 1)  Grid % cells_c(2,c) =-1
          if(j .eq. cj) Grid % cells_c(4,c) =-1
          if(k .eq. 1)  Grid % cells_c(1,c) =-1
          if(k .eq. ck) Grid % cells_c(6,c) =-1

        end do
      end do
    end do

    ! Old number of nodes and cells
    Dom % blocks(b) % n_nodes = Grid % n_nodes
    Dom % blocks(b) % n_cells = Grid % n_cells
    Grid % n_nodes = Grid % n_nodes + ni*nj*nk
    Grid % n_cells = Grid % n_cells + ci*cj*ck

  end do   ! through Dom % blocks

  end subroutine
