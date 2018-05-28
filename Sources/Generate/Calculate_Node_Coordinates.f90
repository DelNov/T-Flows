!==============================================================================!
  subroutine Calculate_Node_Coordinates(dom, grid)
!------------------------------------------------------------------------------!
!   Calculate node coordinates inside the domain, block by block.              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use gen_mod
  use Domain_Mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: grid
!----------------------------------[Calling]-----------------------------------!
  integer :: Is_Line_in_Block
!---------------------------------[Interface]----------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer :: fc, b, bl, i, j, k, n, c, ig
  integer :: l, l1, l2
  integer :: is, js, ks, ie, je, ke
  integer :: ni, nj, nk, ci, cj, ck
  integer :: trans(3,2)
!==============================================================================!

  ! Why on Earth is this?
  do n = 1, grid % max_n_nodes
    grid % xn(n) = HUGE
  end do 

  !------------------------------------!
  !   Calculate the node coordinates   !
  !------------------------------------!
  grid % n_nodes = 0  ! initialize n.o.n.
  grid % n_cells = 0  ! initialize n.o.v.

  do b = 1, size(dom % blocks)

    print '(a38,i7)', '# Generating block:                  ', b
    ni=dom % blocks(b) % resolutions(1)
    nj=dom % blocks(b) % resolutions(2)
    nk=dom % blocks(b) % resolutions(3)   

    ! ( 1 )
    n = grid % n_nodes + ( 1-1)*ni*nj + ( 1-1)*ni +  1
    grid % xn(n) = dom % points(dom % blocks(b) % corners(1)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(1)) % y 
    grid % zn(n) = dom % points(dom % blocks(b) % corners(1)) % z 

    ! ( 2 )
    n = grid % n_nodes + ( 1-1)*ni*nj + ( 1-1)*ni + ni
    grid % xn(n) = dom % points(dom % blocks(b) % corners(2)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(2)) % y
    grid % zn(n) = dom % points(dom % blocks(b) % corners(2)) % z

    ! ( 3 )
    n = grid % n_nodes + ( 1-1)*ni*nj + (nj-1)*ni +  1
    grid % xn(n) = dom % points(dom % blocks(b) % corners(3)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(3)) % y
    grid % zn(n) = dom % points(dom % blocks(b) % corners(3)) % z

    ! ( 4 )
    n = grid % n_nodes + ( 1-1)*ni*nj + (nj-1)*ni + ni
    grid % xn(n) = dom % points(dom % blocks(b) % corners(4)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(4)) % y
    grid % zn(n) = dom % points(dom % blocks(b) % corners(4)) % z

    ! ( 5 ) !
    n = grid % n_nodes + (nk-1)*ni*nj + ( 1-1)*ni +  1
    grid % xn(n) = dom % points(dom % blocks(b) % corners(5)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(5)) % y
    grid % zn(n) = dom % points(dom % blocks(b) % corners(5)) % z

    ! ( 6 ) !
    n = grid % n_nodes + (nk-1)*ni*nj + ( 1-1)*ni + ni
    grid % xn(n) = dom % points(dom % blocks(b) % corners(6)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(6)) % y
    grid % zn(n) = dom % points(dom % blocks(b) % corners(6)) % z

    ! ( 7 ) !
    n = grid % n_nodes + (nk-1)*ni*nj + (nj-1)*ni +  1
    grid % xn(n) = dom % points(dom % blocks(b) % corners(7)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(7)) % y
    grid % zn(n) = dom % points(dom % blocks(b) % corners(7)) % z

    ! ( 8 ) !
    n = grid % n_nodes + (nk-1)*ni*nj + (nj-1)*ni + ni
    grid % xn(n) = dom % points(dom % blocks(b) % corners(8)) % x
    grid % yn(n) = dom % points(dom % blocks(b) % corners(8)) % y
    grid % zn(n) = dom % points(dom % blocks(b) % corners(8)) % z

    !------------------------------!
    !   First on the dom % lines   !
    !    defined point by point    !
    !------------------------------!
    do l=1, size(dom % lines)

      bl = Is_Line_in_Block(dom,                         &
                            dom % lines(l) % points(1),  &
                            dom % lines(l) % points(2),  &
                            b)

      if(bl .eq. b) then

        do n = 1, 8
          if(      dom % lines(l)  % points(1)  &
              .eq. dom % blocks(b) % corners(n) ) l1=n
          if(      dom % lines(l)  % points(2)  &
              .eq. dom % blocks(b) % corners(n) ) l2=n
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
            trans(2,1) =dom % blocks(b) % resolutions(2)
          endif
          if( (l1 .eq. 1).or.(l1 .eq. 3) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =dom % blocks(b) % resolutions(3)
          endif

        ! Line is defined in the -i direction
        else if( (l2-l1) .eq. -1 ) then
          trans(1,1) =dom % blocks(b) % resolutions(1)+1   ! ni from block + 1
          trans(1,2) =-1
          trans(2,2) = 0
          trans(3,2) = 0
          if( (l1 .eq. 2).or.(l1 .eq. 6) ) then
            trans(2,1) =1
          else                         
            trans(2,1) =dom % blocks(b) % resolutions(2)
          endif
          if( (l1 .eq. 2).or.(l1 .eq. 4) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =dom % blocks(b) % resolutions(3)
          endif

        ! Line is defined in the +j direction
        else if( (l2-l1) .eq. +2 ) then
          trans(2,1) = 0
          trans(2,2) =+1
          trans(1,2) = 0
          trans(3,2) = 0 
          if( (l1 .eq. 1).or.(l1 .eq. 5) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =dom % blocks(b) % resolutions(1)
          endif
          if( (l1 .eq. 1).or.(l1 .eq. 2) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =dom % blocks(b) % resolutions(3)
          endif

        ! Line is defined in the -j direction
        else if( (l2-l1) .eq. -2 ) then
          trans(2,1) =dom % blocks(b) % resolutions(2)+1   ! nj from block + 1
          trans(2,2) =-1
          trans(1,2) = 0
          trans(3,2) = 0 
          if( (l1 .eq. 3).or.(l1 .eq. 7) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =dom % blocks(b) % resolutions(1)
          endif
          if( (l1 .eq. 3).or.(l1 .eq. 4) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =dom % blocks(b) % resolutions(3)
          endif

        ! Line is defined in the +k direction
        else if( (l2-l1) .eq. +4 ) then
          trans(3,1) = 0
          trans(3,2) =+1
          trans(1,2) = 0
          trans(2,2) = 0 
          if( (l1 .eq. 1).or.(l1 .eq. 3) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =dom % blocks(b) % resolutions(1)
          end if
          if( (l1 .eq. 1).or.(l1 .eq. 2) ) then
            trans(2,1) =1
          else                         
            trans(2,1) =dom % blocks(b) % resolutions(2)
          endif

        ! Line is defined in the -k direction
        else if( (l2-l1) .eq. -4 ) then
          trans(3,1) =dom % blocks(b) % resolutions(3) + 1  ! nk from block + 1
          trans(3,2) =-1
          trans(1,2) = 0
          trans(2,2) = 0 
          if( (l1 .eq. 5).or.(l1 .eq. 7) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =dom % blocks(b) % resolutions(1)
          endif
          if( (l1 .eq. 5).or.(l1 .eq. 6) ) then
            trans(2,1) =1
          else                         
            trans(2,1) =dom % blocks(b) % resolutions(2)
          endif

        endif ! l1-l2

        ! Line is defined point by point
        if( Approx( dom % lines(l) % weight, 0.0) ) then
          print *, '# Line: ', l
          print *, '# l1= ', l1
          print *, '# l2= ', l2
          do ig=1,dom % lines(l) % resolution
            i=trans(1,1)+trans(1,2)*ig
            j=trans(2,1)+trans(2,2)*ig
            k=trans(3,1)+trans(3,2)*ig

            n = grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i
            grid % xn(n) = dom % lines(l) % x(ig)
            grid % yn(n) = dom % lines(l) % y(ig)
            grid % zn(n) = dom % lines(l) % z(ig)
          end do

        ! Line is defined with a weight factor
        else
          is=trans(1,1)+trans(1,2)
          js=trans(2,1)+trans(2,2)
          ks=trans(3,1)+trans(3,2)
          ie=trans(1,1)+trans(1,2)*dom % lines(l) % resolution
          je=trans(2,1)+trans(2,2)*dom % lines(l) % resolution
          ke=trans(3,1)+trans(3,2)*dom % lines(l) % resolution
          call Distribute_Nodes(dom, grid,                   &
                                b, dom % lines(l) % weight,  &
                                is, js, ks, ie, je, ke)
        endif  

      endif ! if the block contains

    end do ! for the dom % lines

    !-----------!
    !   Lines   !
    !-----------!
    do k=1,nk,nk-1
      do j=1,nj,nj-1
        call Distribute_Nodes(dom, grid,  &
                              b, dom % blocks(b) % weights(1), 1,j,k,ni,j,k)
      end do
    end do

    do k=1,nk,nk-1
      do i=1,ni,ni-1
        call Distribute_Nodes(dom, grid,  &
                              b, dom % blocks(b) % weights(2), i,1,k,i,nj,k)
      end do
    end do

    do j=1,nj,nj-1
      do i=1,ni,ni-1
        call Distribute_Nodes(dom, grid,  &
                              b, dom % blocks(b) % weights(3), i,j,1,i,j,nk)
      end do
    end do

    !------------------------------------------------------------!
    !   Surfaces...                                              !
    !                                                            !
    !   I think this is the propper way to calculate surfaces:   !
    !   it spans the dom % lines in the direction of higher weigh      !
    !------------------------------------------------------------!

    ! I (k=1) 
    fc = 1   ! face index
    k = 1
    if( .not. Approx(dom % blocks(b) % face_weights(fc,1),1.0 ) ) then
      do j=1,nj
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    else ! dom % lines in the j direction
      do i=1,ni
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    endif

    ! VI (k=nk)
    fc = 6   ! face index
    k = nk
    if( .not. Approx(dom % blocks(b) % face_weights(fc,1),1.0 ) ) then
     do j=1,nj
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    else ! dom % lines in the j direction
      do i=1,ni
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    endif

    ! V (i=1)
    fc = 5   ! face index
    i = 1
    if( .not. Approx(dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do j=1,nj
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,3),  &
                              i,j,1,i,j,nk)
      end do
    else ! dom % lines in the j direction
      do k=1,nk
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    end if 

    ! III (i=ni)
    fc = 3   ! face index
    i = ni
    if( .not. Approx(dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do j=1,nj
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,3),  & 
                              i,j,1,i,j,nk)
      end do
    else ! dom % lines in the j direction
      do k=1,nk
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,2),  &
                              i,1,k,i,nj,k)
      end do
    end if 

    ! II (j=1)       
    fc = 2   ! face index
    j = 1
    if( .not. Approx(dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do i=1,ni
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,3),  &
                              i,j,1,i,j,nk)
      end do
    else ! dom % lines in the i direction
      do k=1,nk
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    endif

    ! IV (j=nj)       
    fc = 4   ! face index
    j = nj
    if( .not. Approx(dom % blocks(b) % face_weights(fc,3),1.0 ) ) then
      do i=1,ni
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,3),  &
                              i,j,1,i,j,nk)
      end do
    else ! dom % lines in the i direction
      do k=1,nk
        call Distribute_Nodes(dom, grid,                                &
                              b, dom % blocks(b) % face_weights(fc,1),  &
                              1,j,k,ni,j,k)
      end do
    endif

    !-------------!
    !   Volumes   !
    !-------------!
    if( .not. Approx( dom % blocks(b) % weights(3), 1.0 ) ) then
      do i=1,ni
        do j=1,nj
          call Distribute_Nodes(dom, grid,  &
                                b, dom % blocks(b) % weights(3), i,j,1,i,j,nk)
        end do
      end do
    else if( .not. Approx( dom % blocks(b) % weights(1), 1.0 ) ) then
      do k=1,nk
        do j=1,nj
          call Distribute_Nodes(dom, grid,  &
                                b, dom % blocks(b) % weights(1), 1,j,k,ni,j,k)
        end do
      end do
    else if( .not. Approx( dom % blocks(b) % weights(2), 1.0 ) ) then
      do k=1,nk
        do i=1,ni
          call Distribute_Nodes(dom, grid,  &
                                b, dom % blocks(b) % weights(2), i,1,k,i,nj,k)
        end do
      end do
    else

      do i=1,ni
        do j=1,nj
          do k=1,nk
            n = grid % n_nodes+(k-1)*ni*nj + (j-1)*ni + i
            call Laplace(dom, grid, b, i, j, k,            &
                         ONE_THIRD, ONE_THIRD, ONE_THIRD,  &
                         ONE_THIRD, ONE_THIRD, ONE_THIRD,  &
                         ONE_THIRD, ONE_THIRD, ONE_THIRD)
          end do
        end do
      end do
    endif

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
          c = grid % n_cells + (k-1)*ci*cj + (j-1)*ci + i ! cell 
          n = grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i ! 1st node

          ! Nodes
          grid % cells_n(1,c) = n
          grid % cells_n(2,c) = n+1
          grid % cells_n(3,c) = n+ni
          grid % cells_n(4,c) = n+ni+1
          grid % cells_n(5,c) = grid % cells_n(1,c)+ni*nj
          grid % cells_n(6,c) = grid % cells_n(2,c)+ni*nj
          grid % cells_n(7,c) = grid % cells_n(3,c)+ni*nj
          grid % cells_n(8,c) = grid % cells_n(4,c)+ni*nj

          ! Neighbours
          grid % cells_c(1,c) = c-ci*cj
          grid % cells_c(2,c) = c-ci 
          grid % cells_c(3,c) = c+1
          grid % cells_c(4,c) = c+ci
          grid % cells_c(5,c) = c-1
          grid % cells_c(6,c) = c+ci*cj

          ! This value (-1) is also the default boundary marker
          if(i .eq. 1)  grid % cells_c(5,c) =-1
          if(i .eq. ci) grid % cells_c(3,c) =-1
          if(j .eq. 1)  grid % cells_c(2,c) =-1
          if(j .eq. cj) grid % cells_c(4,c) =-1
          if(k .eq. 1)  grid % cells_c(1,c) =-1
          if(k .eq. ck) grid % cells_c(6,c) =-1

        end do
      end do
    end do

    ! Old number of nodes and cells
    dom % blocks(b) % n_nodes = grid % n_nodes  
    dom % blocks(b) % n_cells = grid % n_cells  
    grid % n_nodes = grid % n_nodes + ni*nj*nk
    grid % n_cells = grid % n_cells + ci*cj*ck

  end do   ! through dom % blocks 

  end subroutine
