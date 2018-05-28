!==============================================================================!
  subroutine Connect_Copy(dom, grid)
!------------------------------------------------------------------------------!
!   Solve the cell connectivity for copy boundary conditions.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod
  use Domain_Mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n, p, n_cop                ! counters
  integer :: b1, b2                           ! block 1 and 2
  integer :: f1, f2                           ! faces of block 1 and 2
  integer :: n11,n12,n13,n14,n21,n22,n23,n24  ! global node numbers
  integer :: p11,p12,p13,p14,p21,p22,p23,p24  ! global node numbers
  integer :: l11,l12,l13,l14,l21,l22,l23,l24  ! local  node numbers
  integer :: i1, j1, i2, j2, k1, k2           ! directions in dom % blocks
  integer :: ig, jg, nig, njg                 ! generic plane 
  integer :: ci1, cj1, ck1, ci2, cj2, ck2     ! resolution of dom % blocks
  integer :: c1, c2                           ! cells from block 1, 2
  integer :: trans1(3,3), trans2(3,3)
!==============================================================================!

  !-----------------------------------------------------!
  !   Search through all block and all of their faces   !
  !-----------------------------------------------------!
  n_cop = 0       

  do p = 1, n_copy_cond    
    do b2 = 1, size(dom % blocks)
      do b1 = 1, size(dom % blocks)
        do f2 = 1, 6    ! faces of the second block
          do f1 = 1, 6  ! faces of the first block

            ! Initialize the transformation matrixes             
            do i=1,3
              do j=1,3
                trans1(i,j)=0
                trans2(i,j)=0
              end do
            end do

            n11 = dom % blocks(b1) % faces(f1, 1)
            n12 = dom % blocks(b1) % faces(f1, 2)
            n13 = dom % blocks(b1) % faces(f1, 3)
            n14 = dom % blocks(b1) % faces(f1, 4) 
            n21 = dom % blocks(b2) % faces(f2, 1)
            n22 = dom % blocks(b2) % faces(f2, 2)
            n23 = dom % blocks(b2) % faces(f2, 3)
            n24 = dom % blocks(b2) % faces(f2, 4)

            p11=copy_cond(p, 1)
            p12=copy_cond(p, 2)
            p13=copy_cond(p, 3)
            p14=copy_cond(p, 4) 
            p21=copy_cond(p, 5)
            p22=copy_cond(p, 6)
            p23=copy_cond(p, 7)
            p24=copy_cond(p, 8)

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
                if(dom % blocks(b1) % corners(n) .eq. p11) l11 = n
                if(dom % blocks(b1) % corners(n) .eq. p12) l12 = n
                if(dom % blocks(b1) % corners(n) .eq. p13) l13 = n
                if(dom % blocks(b1) % corners(n) .eq. p14) l14 = n
                if(dom % blocks(b2) % corners(n) .eq. p21) l21 = n
                if(dom % blocks(b2) % corners(n) .eq. p22) l22 = n
                if(dom % blocks(b2) % corners(n) .eq. p23) l23 = n
                if(dom % blocks(b2) % corners(n) .eq. p24) l24 = n
              end do

              write(*, *) '# Copy_cond dom % blocks: ', b1, b2

              ! Direction ig, block 1
              if((l14-l11) .eq. +1) then
                nig = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,2)=+1
              elseif((l14-l11) .eq. +2) then
                nig = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,2)=+1
              elseif((l14-l11) .eq. +4) then 
                nig = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,2)=+1
              elseif((l14-l11) .eq. -1) then 
                nig = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,1)=nig
                trans1(1,2)=-1
              elseif((l14-l11) .eq. -2) then 
                nig = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,1)=nig
                trans1(2,2)=-1
              elseif((l14-l11) .eq. -4) then 
                nig = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,1)=nig
                trans1(3,2)=-1
              endif

              ! Direction jg, block 1 
              if((l12-l11) .eq. +1) then 
                njg = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,3)=+1
              elseif((l12-l11) .eq. +2) then
                njg = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,3)=+1
              elseif((l12-l11) .eq. +4) then
                njg = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,3)=+1
              elseif((l12-l11) .eq. -1) then
                njg = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,1)=njg
                trans1(1,3)=-1
              elseif((l12-l11) .eq. -2) then
                njg = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,1)=njg
                trans1(2,3)=-1
              elseif((l12-l11) .eq. -4) then
                njg = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,1)=njg
                trans1(3,3)=-1
              endif

              ! Direction ig, block 2
              if((l24-l21) .eq. +1) then
                nig = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,2)=+1
              elseif((l24-l21) .eq. +2) then
                nig = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,2)=+1
              elseif((l24-l21) .eq. +4) then 
                nig = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,2)=+1
              elseif((l24-l21) .eq. -1) then 
                nig = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,1)=nig
                trans2(1,2)=-1
              elseif((l24-l21) .eq. -2) then 
                nig = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,1)=nig
                trans2(2,2)=-1
              elseif((l24-l21) .eq. -4) then 
                nig = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,1)=nig
                trans2(3,2)=-1
              endif

              ! Direction jg, block 2 
              if((l22-l21) .eq. +1) then 
                njg = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,3)=+1
              elseif((l22-l21) .eq. +2) then
                njg = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,3)=+1
              elseif((l22-l21) .eq. +4) then
                njg = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,3)=+1
              elseif((l22-l21) .eq. -1) then
                njg = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,1)=njg
                trans2(1,3)=-1
              elseif((l22-l21) .eq. -2) then
                njg = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,1)=njg
                trans2(2,3)=-1
              elseif((l22-l21) .eq. -4) then
                njg = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,1)=njg
                trans2(3,3)=-1
              endif

              ! Set the constant directions
              if(f1 .eq. 1) trans1(3,1)=1
              if(f1 .eq. 2) trans1(2,1)=1
              if(f1 .eq. 3) trans1(1,1)=dom % blocks(b1) % resolutions(1)-1
              if(f1 .eq. 4) trans1(2,1)=dom % blocks(b1) % resolutions(2)-1
              if(f1 .eq. 5) trans1(1,1)=1
              if(f1 .eq. 6) trans1(3,1)=dom % blocks(b1) % resolutions(3)-1

              if(f2 .eq. 1) trans2(3,1)=1
              if(f2 .eq. 2) trans2(2,1)=1
              if(f2 .eq. 3) trans2(1,1)=dom % blocks(b2) % resolutions(1)-1
              if(f2 .eq. 4) trans2(2,1)=dom % blocks(b2) % resolutions(2)-1
              if(f2 .eq. 5) trans2(1,1)=1
              if(f2 .eq. 6) trans2(3,1)=dom % blocks(b2) % resolutions(3)-1 

              ! Finally conect the two copy boundaries
              ! now, connect the cells in the interior
              do jg=1,njg-1              ! through volumes only
                do ig=1,nig-1            ! through volumes only
                  n_cop = n_cop + 1
                  ci1=dom % blocks(b1) % resolutions(1)-1
                  cj1=dom % blocks(b1) % resolutions(2)-1
                  ck1=dom % blocks(b1) % resolutions(3)-1
                  ci2=dom % blocks(b2) % resolutions(1)-1
                  cj2=dom % blocks(b2) % resolutions(2)-1
                  ck2=dom % blocks(b2) % resolutions(3)-1
                  i1 = trans1(1,1)+trans1(1,2)*ig+trans1(1,3)*jg
                  j1 = trans1(2,1)+trans1(2,2)*ig+trans1(2,3)*jg
                  k1 = trans1(3,1)+trans1(3,2)*ig+trans1(3,3)*jg 
                  i2 = trans2(1,1)+trans2(1,2)*ig+trans2(1,3)*jg
                  j2 = trans2(2,1)+trans2(2,2)*ig+trans2(2,3)*jg
                  k2 = trans2(3,1)+trans2(3,2)*ig+trans2(3,3)*jg
                  c1 = dom % blocks(b1) % n_cells  &
                     + (k1-1)*ci1*cj1 + (j1-1)*ci1 + i1
                  c2 = dom % blocks(b2) % n_cells  &
                     + (k2-1)*ci2*cj2 + (j2-1)*ci2 + i2
                  grid % bnd_cond % copy_c(c2) = c1 ! allway copy from c1 to c2
                  grid % bnd_cond % copy_s(1, n_cop) = c1
                  grid % bnd_cond % copy_s(2, n_cop) = c2
                end do
              end do        

            end if  ! are they connected ? 

          end do    ! f1
        end do      ! f2
      end do        ! b1
    end do          ! b2 
  end do            ! p copies 

  n_copy_cond = n_cop

  end subroutine
