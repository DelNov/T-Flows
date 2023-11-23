!==============================================================================!
  subroutine Connect_Blocks(Dom, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to establish cell connectivity across different
!>  blocks in a grid after they have been generated.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initialization: Sets up a mapping array (Grid % new_n) to keep track of  !
!     new node indices after connecting blocks and initializes the number of   !
!     deleted nodes (del).                                                     !
!   * Early exit for single block: If there's only one block in the domain,    !
!     the subroutine returns early as there's no need for inter-block          !
!     connectivity.                                                            !
!   * Iterating through block faces: It iterates over pairs of blocks,         !
!     checking each face of one block against every face of the other to       !
!     find matching faces that need to be connected.                           !
!   * Node and cell connectivity: Once matching faces are identified, the      !
!     subroutine establishes node-to-node connectivity between these faces.    !
!     It uses transformation matrices (trans1, trans2) to map the coordinates  !
!     from one block to another, ensuring that nodes on the shared face of     !
!     two blocks are correctly aligned and connected.                          !
!   * Updating node and cell information: After connecting the blocks, it      !
!     updates the global node and cell information in the Grid, adjusting      !
!     node numbers to account for merged nodes and ensuring that the cell      !
!     connectivity reflects the new structure.                                 !
!   * Final adjustments and cleanup: The subroutine finalizes the process by   !
!     adjusting the global node count to reflect the new configuration and     !
!     updating the node references in the cell structure (Grid % cells_n)      !
!     to match the new node indices.                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type) :: Dom   !! domain in which the grid is being generated
  type(Grid_Type)    :: Grid  !! grid being generated
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n                         ! counters
  integer :: b1, b2                          ! block 1 and 2
  integer :: f1, f2                          ! faces of block 1 and 2
  integer :: n11,n12,n13,n14,n21,n22,n23,n24 ! global node numbers
  integer :: l11,l12,l13,l14,l21,l22,l23,l24 ! local  node numbers
  integer :: g1, g2, g3, g4                  ! generic points
  integer :: i1, j1, i2, j2, k1, k2          ! directions in Dom % blocks
  integer :: ig, jg, nig, njg                ! generic plane
  integer :: ci1, cj1, ck1, ci2, cj2, ck2    ! resolution of Dom % blocks
  integer :: c1, c2                          ! cells from block 1, 2
  integer :: ni1, nj1, nk1, ni2, nj2, nk2    ! resolution of Dom % blocks
  integer :: n1, n2                          ! from block 1, 2
  integer :: trans1(3,3), trans2(3,3)
  integer :: del                             ! number of deleted nodes
!==============================================================================!

  ! Initialize the new_n array
  do n = 1, Grid % n_nodes
    Grid % new_n(n) = n
  end do

  ! If number of Dom % blocks is equal to one, there is nothing to do
  if(size(Dom % blocks) .eq. 1) return

  ! Initialize the number of deleted nodes
  del=0

  ni1=0; nj1=0; nk1=0; ni2=0; nj2=0; nk2=0

  !-----------------------------------------------------!
  !   Search through all block and all of their faces   !
  !-----------------------------------------------------!
  do b2 = 2, size(Dom % blocks)
    do b1 = 1, b2-1
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

          ! Check if they are connected
          if( ((n11 .eq. n21).and.(n13 .eq. n23)) .or.  &
              ((n11 .eq. n24).and.(n13 .eq. n22)) .or.  &
              ((n11 .eq. n23).and.(n13 .eq. n21)) .or.  &
              ((n11 .eq. n22).and.(n13 .eq. n24)) ) then

            ! Define generic surface (g1-g4 are in essence not needed)
            g1=n11
            g2=n12
            g3=n13
            g4=n14

            ! Find local nodes (1-8) from blocks 1 and 2 on generic surface
            do n = 1, 8
              if(Dom % blocks(b1) % corners(n) .eq. g1) l11=n
              if(Dom % blocks(b2) % corners(n) .eq. g1) l21=n
              if(Dom % blocks(b1) % corners(n) .eq. g2) l12=n
              if(Dom % blocks(b2) % corners(n) .eq. g2) l22=n
              if(Dom % blocks(b1) % corners(n) .eq. g3) l13=n
              if(Dom % blocks(b2) % corners(n) .eq. g3) l23=n
              if(Dom % blocks(b1) % corners(n) .eq. g4) l14=n
              if(Dom % blocks(b2) % corners(n) .eq. g4) l24=n
            end do

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

            ! Finally conect the two Dom % blocks
            do jg=1,njg-1              ! through cells only
              do ig=1,nig-1            ! through cells only
                ci1=Dom % blocks(b1) % resolutions(1)-1
                cj1=Dom % blocks(b1) % resolutions(2)-1
                ck1=Dom % blocks(b1) % resolutions(3)-1
                ci2=Dom % blocks(b2) % resolutions(1)-1
                cj2=Dom % blocks(b2) % resolutions(2)-1
                ck2=Dom % blocks(b2) % resolutions(3)-1
                i1 = trans1(1,1) + trans1(1,2)*ig + trans1(1,3)*jg
                j1 = trans1(2,1) + trans1(2,2)*ig + trans1(2,3)*jg
                k1 = trans1(3,1) + trans1(3,2)*ig + trans1(3,3)*jg
                i2 = trans2(1,1) + trans2(1,2)*ig + trans2(1,3)*jg
                j2 = trans2(2,1) + trans2(2,2)*ig + trans2(2,3)*jg
                k2 = trans2(3,1) + trans2(3,2)*ig + trans2(3,3)*jg
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
                Grid % new_n(n2) = Grid % new_n(n1)
              end do
            end do

          end if  ! are they connected ?

        end do    ! f1
      end do      ! f2
    end do        ! b1

    ! Update node numbers
    do n = Dom % blocks(b2) % n_nodes + 1,  &
           Dom % blocks(b2) % n_nodes + ni2*nj2*nk2
      if(Grid % new_n(n) .ne. n) del = del + 1
      if(Grid % new_n(n) .eq. n) Grid % new_n(n) = Grid % new_n(n) - del
    end do

  end do          ! b2

  do n = 1, Grid % n_nodes
    Grid % xn(Grid % new_n(n)) = Grid % xn(n)
    Grid % yn(Grid % new_n(n)) = Grid % yn(n)
    Grid % zn(Grid % new_n(n)) = Grid % zn(n)
  end do

  Grid % n_nodes = Grid % n_nodes - del

  ! Skip the merged points in the node() structure
  do i = 1, Grid % n_cells
    do n = 1, 8
      Grid % cells_n(n,i) = Grid % new_n(Grid % cells_n(n, i))
    end do
  end do

  end subroutine
