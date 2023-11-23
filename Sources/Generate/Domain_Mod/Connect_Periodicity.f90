!==============================================================================!
  subroutine Connect_Periodicity(Dom, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to establish cell connectivity for periodic
!>  boundary conditions in a computational grid.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initialization of twin nodes: Initializes a mapping array (twin_n)       !
!     that tracks nodes across periodic boundaries.                            !
!   * Iterating through blocks and faces: Loops through all blocks and their   !
!     faces to identify potential periodic connections.                        !
!   * Transformation matrices setup: Initializes transformation matrices       !
!     (trans1, trans2) to map node coordinates between periodic faces.         !
!   * Checking for periodic connectivity: Compares nodes on faces of different !
!     blocks to determine if they are periodic counterparts, based on          !
!     predefined periodic conditions.                                          !
!   * Establishing cell connectivity: For faces that are identified as         !
!     periodic counterparts, it sets up cell connectivity across the periodic  !
!     boundary, mapping cells from one block to the corresponding cells in the !
!     periodic partner block.                                                  !
!   * Node connectivity across periodic boundaries: Establishes node-to-node   !
!     connections across the periodic boundaries, ensuring that each node on   !
!     a periodic face is correctly paired with its counterpart.                !
!   * Twin of my twin logic: Implements a logic to ensure that if a node is a  !
!     twin of another node, which in turn is a twin of a third node, then the  !
!     first node is also considered a twin of the third node. This step        !
!     ensures comprehensive connectivity across all periodic nodes.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type) :: Dom   !! domain in which the grid is being generated
  type(Grid_Type)    :: Grid  !! grid being generated
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n, p                      ! counters
  integer :: b1, b2                          ! block 1 and 2
  integer :: f1, f2                          ! faces of block 1 and 2
  integer :: n11,n12,n13,n14,n21,n22,n23,n24 ! global node numbers
  integer :: p11,p12,p13,p14,p21,p22,p23,p24 ! global node numbers
  integer :: l11,l12,l13,l14,l21,l22,l23,l24 ! local  node numbers
  integer :: i1, j1, i2, j2, k1, k2          ! directions in Dom % blocks
  integer :: ig, jg, nig, njg                ! generic plane
  integer :: ci1, cj1, ck1, ci2, cj2, ck2    ! resolution of Dom % blocks
  integer :: c1, c2                          ! cells from block 1, 2
  integer :: ni1, nj1, nk1, ni2, nj2, nk2    ! resolution of Dom % blocks
  integer :: n1, n2                          !  from block 1, 2
  integer :: n3, i3, new
  integer :: trans1(3,3), trans2(3,3)
!==============================================================================!

  ! Initialise twin_n
  do n = 1, Grid % max_n_nodes
    twin_n(n,0) = 0
  end do

  !-----------------------------------------------------!
  !   Search through all block and all of their faces   !
  !-----------------------------------------------------!
  do p = 1, n_periodic_cond
    do b2 = 1, size(Dom % blocks)
      do b1 = 1, size(Dom % blocks)
        do f2 = 1, 6    ! faces of the second block
          do f1 = 1, 6  ! faces of the first block

            ! Initialize the transformation matrixes
            do i=1,3
              do j=1,3
                trans1(i,j)=0
                trans2(i,j)=0
              end do
            end do

            n11 = Dom % blocks(b1) % faces(f1, 1)
            n12 = Dom % blocks(b1) % faces(f1, 2)
            n13 = Dom % blocks(b1) % faces(f1, 3)
            n14 = Dom % blocks(b1) % faces(f1, 4)
            n21 = Dom % blocks(b2) % faces(f2, 1)
            n22 = Dom % blocks(b2) % faces(f2, 2)
            n23 = Dom % blocks(b2) % faces(f2, 3)
            n24 = Dom % blocks(b2) % faces(f2, 4)

            p11=periodic_cond(p, 1)
            p12=periodic_cond(p, 2)
            p13=periodic_cond(p, 3)
            p14=periodic_cond(p, 4)
            p21=periodic_cond(p, 5)
            p22=periodic_cond(p, 6)
            p23=periodic_cond(p, 7)
            p24=periodic_cond(p, 8)

            ! Check if they are connected
            if( ( ((n11 .eq. p11).and.(n13 .eq. p13)) .or.                &
                  ((n11 .eq. p14).and.(n13 .eq. p12)) .or.                &
                  ((n11 .eq. p13).and.(n13 .eq. p11)) .or.                &
                  ((n11 .eq. p12).and.(n13 .eq. p14)) )                   &
                                 .and.                                    &
                ( ((n21 .eq. p21).and.(n23 .eq. p23)) .or.                &
                  ((n21 .eq. p24).and.(n23 .eq. p22)) .or.                &
                  ((n21 .eq. p23).and.(n23 .eq. p21)) .or.                &
                  ((n21 .eq. p22).and.(n23 .eq. p24)) ) ) then

              ! Find local nodes (1-8) from blocks 1 and 2 on generic surface
              do n=1,8
                if(Dom % blocks(b1) % corners(n) .eq. p11) l11=n
                if(Dom % blocks(b1) % corners(n) .eq. p12) l12=n
                if(Dom % blocks(b1) % corners(n) .eq. p13) l13=n
                if(Dom % blocks(b1) % corners(n) .eq. p14) l14=n
                if(Dom % blocks(b2) % corners(n) .eq. p21) l21=n
                if(Dom % blocks(b2) % corners(n) .eq. p22) l22=n
                if(Dom % blocks(b2) % corners(n) .eq. p23) l23=n
                if(Dom % blocks(b2) % corners(n) .eq. p24) l24=n
              end do

               print '(a31,2i7)', '# Periodicity between blocks: ', b1, b2

              ! Direction ig, block 1
              if((l14-l11) .eq. +1) then
                nig = Dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,2)=+1
              elseif((l14-l11) .eq. +2) then
                nig = Dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,2)=+1
              elseif((l14-l11) .eq. +4) then
                nig = Dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,2)=+1
              elseif((l14-l11) .eq. -1) then
                nig = Dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,1)=nig
                trans1(1,2)=-1
              elseif((l14-l11) .eq. -2) then
                nig = Dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,1)=nig
                trans1(2,2)=-1
              elseif((l14-l11) .eq. -4) then
                nig = Dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,1)=nig
                trans1(3,2)=-1
              end if

              ! Direction jg, block 1
              if((l12-l11) .eq. +1) then
                njg = Dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,3)=+1
              elseif((l12-l11) .eq. +2) then
                njg = Dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,3)=+1
              elseif((l12-l11) .eq. +4) then
                njg = Dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,3)=+1
              elseif((l12-l11) .eq. -1) then
                njg = Dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,1)=njg
                trans1(1,3)=-1
              elseif((l12-l11) .eq. -2) then
                njg = Dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,1)=njg
                trans1(2,3)=-1
              elseif((l12-l11) .eq. -4) then
                njg = Dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,1)=njg
                trans1(3,3)=-1
              end if

              ! Direction ig, block 2
              if((l24-l21) .eq. +1) then
                nig = Dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,2)=+1
              elseif((l24-l21) .eq. +2) then
                nig = Dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,2)=+1
              elseif((l24-l21) .eq. +4) then
                nig = Dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,2)=+1
              elseif((l24-l21) .eq. -1) then
                nig = Dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,1)=nig
                trans2(1,2)=-1
              elseif((l24-l21) .eq. -2) then
                nig = Dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,1)=nig
                trans2(2,2)=-1
              elseif((l24-l21) .eq. -4) then
                nig = Dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,1)=nig
                trans2(3,2)=-1
              end if

              ! Direction jg, block 2
              if((l22-l21) .eq. +1) then
                njg = Dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,3)=+1
              elseif((l22-l21) .eq. +2) then
                njg = Dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,3)=+1
              elseif((l22-l21) .eq. +4) then
                njg = Dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,3)=+1
              elseif((l22-l21) .eq. -1) then
                njg = Dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,1)=njg
                trans2(1,3)=-1
              elseif((l22-l21) .eq. -2) then
                njg = Dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,1)=njg
                trans2(2,3)=-1
              elseif((l22-l21) .eq. -4) then
                njg = Dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,1)=njg
                trans2(3,3)=-1
              end if

              ! Set the constant directions
              if(f1 .eq. 1) trans1(3,1)=1
              if(f1 .eq. 2) trans1(2,1)=1
              if(f1 .eq. 3) trans1(1,1)=Dom % blocks(b1) % resolutions(1)-1
              if(f1 .eq. 4) trans1(2,1)=Dom % blocks(b1) % resolutions(2)-1
              if(f1 .eq. 5) trans1(1,1)=1
              if(f1 .eq. 6) trans1(3,1)=Dom % blocks(b1) % resolutions(3)-1

              if(f2 .eq. 1) trans2(3,1)=1
              if(f2 .eq. 2) trans2(2,1)=1
              if(f2 .eq. 3) trans2(1,1)=Dom % blocks(b2) % resolutions(1)-1
              if(f2 .eq. 4) trans2(2,1)=Dom % blocks(b2) % resolutions(2)-1
              if(f2 .eq. 5) trans2(1,1)=1
              if(f2 .eq. 6) trans2(3,1)=Dom % blocks(b2) % resolutions(3)-1

              ! Finally conect the two periodic boundaries
              do jg=1,njg-1              ! through cells only
                do ig=1,nig-1            ! through cells only
                  ci1=Dom % blocks(b1) % resolutions(1)-1
                  cj1=Dom % blocks(b1) % resolutions(2)-1
                  ck1=Dom % blocks(b1) % resolutions(3)-1
                  ci2=Dom % blocks(b2) % resolutions(1)-1
                  cj2=Dom % blocks(b2) % resolutions(2)-1
                  ck2=Dom % blocks(b2) % resolutions(3)-1
                  i1 = trans1(1,1)+trans1(1,2)*ig+trans1(1,3)*jg
                  j1 = trans1(2,1)+trans1(2,2)*ig+trans1(2,3)*jg
                  k1 = trans1(3,1)+trans1(3,2)*ig+trans1(3,3)*jg
                  i2 = trans2(1,1)+trans2(1,2)*ig+trans2(1,3)*jg
                  j2 = trans2(2,1)+trans2(2,2)*ig+trans2(2,3)*jg
                  k2 = trans2(3,1)+trans2(3,2)*ig+trans2(3,3)*jg
                  c1 = Dom % blocks(b1) % n_cells  &
                     + (k1-1)*ci1*cj1 + (j1-1)*ci1 + i1
                  c2 = Dom % blocks(b2) % n_cells  &
                     + (k2-1)*ci2*cj2 + (j2-1)*ci2 + i2
                  Grid % cells_c(f1,c1) = c2
                  Grid % cells_c(f2,c2) = c1
                end do
              end do

              ! Modify the transformation matrices for nodal connection
              if(trans1(1,1)  > 1) trans1(1,1)=trans1(1,1)+1
              if(trans1(2,1)  > 1) trans1(2,1)=trans1(2,1)+1
              if(trans1(3,1)  > 1) trans1(3,1)=trans1(3,1)+1
              if(trans2(1,1)  > 1) trans2(1,1)=trans2(1,1)+1
              if(trans2(2,1)  > 1) trans2(2,1)=trans2(2,1)+1
              if(trans2(3,1)  > 1) trans2(3,1)=trans2(3,1)+1

              ! Conect the nodes
              do jg=1,njg                ! through nodes
                do ig=1,nig              ! through nodes
                  ni1=Dom % blocks(b1) % resolutions(1)
                  nj1=Dom % blocks(b1) % resolutions(2)
                  nk1=Dom % blocks(b1) % resolutions(3)
                  ni2=Dom % blocks(b2) % resolutions(1)
                  nj2=Dom % blocks(b2) % resolutions(2)
                  nk2=Dom % blocks(b2) % resolutions(3)
                  i1 = trans1(1,1) + trans1(1,2)*ig + trans1(1,3)*jg
                  j1 = trans1(2,1) + trans1(2,2)*ig + trans1(2,3)*jg
                  k1 = trans1(3,1) + trans1(3,2)*ig + trans1(3,3)*jg
                  i2 = trans2(1,1) + trans2(1,2)*ig + trans2(1,3)*jg
                  j2 = trans2(2,1) + trans2(2,2)*ig + trans2(2,3)*jg
                  k2 = trans2(3,1) + trans2(3,2)*ig + trans2(3,3)*jg
                  n1 = Dom % blocks(b1) % n_nodes  &
                     + (k1-1)*ni1*nj1 + (j1-1)*ni1 + i1
                  n2 = Dom % blocks(b2) % n_nodes  &
                     + (k2-1)*ni2*nj2 + (j2-1)*ni2 + i2
                  n1 = Grid % new_n(n1)
                  n2 = Grid % new_n(n2)

                  ! Check if they are already connected
                  do n=1, twin_n(n1,0)
                    if(Gen_Mod_Are_Nodes_Twins(n1,n2)) goto 1
                  end do

                  ! If they were not, connect them
                  twin_n(n1,0)=twin_n(n1,0)+1
                  twin_n(n1,twin_n(n1,0))=n2
                  twin_n(n2,0)=twin_n(n2,0)+1
                  twin_n(n2,twin_n(n2,0))=n1

1               end do  ! jg
              end do    ! ig

            end if  ! are they connected ?

          end do    ! f1
        end do      ! f2
      end do        ! b1
    end do          ! b2
  end do            ! p periods

  !---------------------!
  !   Twin of my twin   !
  !   is also my twin   !
  !---------------------!
  do n1=1,Grid % n_nodes
    do i1=1,twin_n(n1,0)
      n2=twin_n(n1,i1)
      do i2=1,twin_n(n2,0)
        n3=twin_n(n2,i2)   ! twins from n2
        new=n3
        do i3=1,twin_n(n1,0)
          if( (twin_n(n1,i3) .eq. n3) .or. (n3 .eq. n1) ) new=0
        end do
        if(new .eq. n3) then
          twin_n(n1,0)=twin_n(n1,0)+1
          twin_n(n1,twin_n(n1,0))=n3
        end if
      end do
    end do
  end do

  end subroutine
