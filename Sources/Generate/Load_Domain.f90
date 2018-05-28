!==============================================================================!
  subroutine Load_Domain(dom, grid)
!------------------------------------------------------------------------------!
!   Reads: .dom file                                                           !
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use gen_mod
  use Tokenizer_Mod
  use Domain_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: grid
!----------------------------------[Calling]-----------------------------------!
  real :: Tet_Volume   
!-----------------------------------[Locals]-----------------------------------!
  integer           :: b, i, l, s, fc, n, n1,n2,n3,n4
  integer           :: n_faces_check, n_nodes_check
  integer           :: ni, nj, nk, npnt, nsurf
  character(len=12) :: dum
  character(len=80) :: domain_name
  character(len=12) :: answer
  real              :: xt(8), yt(8), zt(8)
  integer           :: face_nodes(6,4)
!==============================================================================!
  data face_nodes / 1, 1, 2, 4, 3, 5,         &
                    2, 5, 6, 8, 7, 7,         &
                    4, 6, 8, 7, 5, 8,         &
                    3, 2, 4, 3, 1, 6  /
!------------------------------------------------------------------------------!

  print *, '#========================================'
  print *, '# Input problem name: (without extension)'
  print *, '#----------------------------------------'
  call Tokenizer_Mod_Read_Line(5) 
  read(line % tokens(1), *) problem_name

  domain_name = problem_name
  domain_name(len_trim(problem_name)+1:len_trim(problem_name)+4) = '.dom'
  write(*, *) '# Now reading the file: ', domain_name
  open(9, file=domain_name)

  !-----------------------------------------------------------------!
  !   Max. number of nodes (cells), boundary faces and cell faces   ! 
  !-----------------------------------------------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) grid % max_n_nodes
  read(line % tokens(2), *) grid % max_n_bnd_cells
  read(line % tokens(3), *) grid % max_n_faces  

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  print *, '# Allocating memory for: ' 
  print *, '#', grid % max_n_nodes,     ' nodes and cells' 
  print *, '#', grid % max_n_bnd_cells, ' boundary cells'         
  print *, '#', grid % max_n_faces,     ' cell faces' 

  allocate (grid % bnd_cond % copy_c(-grid % max_n_bnd_cells:grid % max_n_nodes))
  grid % bnd_cond % copy_c = 0
  allocate (grid % bnd_cond % copy_s(2,grid % max_n_bnd_cells))
  grid % bnd_cond % copy_s = 0    

  allocate (grid % bnd_cond % color(-grid % max_n_bnd_cells-1:-1))
  grid % bnd_cond % color = 0

  ! Variables in Grid_Mod
  call Grid_Mod_Allocate_Nodes(grid,  &
                               grid % max_n_nodes) 

  call Grid_Mod_Allocate_Cells(grid,                    &
                               grid % max_n_nodes,      &
                               grid % max_n_bnd_cells) 

  call Grid_Mod_Allocate_Faces(grid,                &
                               grid % max_n_faces,  &
                               0)                      ! no shadow faces

  ! Variables declared in gen_mod.h90:
  allocate (face_c_to_c(grid % max_n_faces,2))
  face_c_to_c=0 

  ! Variables for renumbering
  allocate (new_n(-grid % max_n_bnd_cells:grid % max_n_nodes))
  new_n=0
  allocate (new_c(-grid % max_n_bnd_cells:grid % max_n_nodes))
  new_c=0
  allocate (new_f( grid % max_n_faces))
  new_f=0
  allocate (cell_marked(-grid % max_n_bnd_cells:grid % max_n_nodes))
  cell_marked = .false.

  allocate (twin_n(grid % max_n_nodes,0:8))
  twin_n=0

  allocate (level(grid % max_n_nodes))
  level=0

  print *, '# Allocation successfull !'

  !-------------!
  !   Corners   !
  !-------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) dom % n_points  ! number of points

  call Domain_Mod_Allocate_Points(dom, dom % n_points)

  do i = 1, dom % n_points
    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(2),*) dom % points(i) % x
    read(line % tokens(3),*) dom % points(i) % y
    read(line % tokens(4),*) dom % points(i) % z
  end do

  !------------!
  !   Blocks   !
  !------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) dom % n_blocks  ! number of blocks 

  call Domain_Mod_Allocate_Blocks(dom, dom % n_blocks)

  ! Initialize weights 
  do b=1, dom % n_blocks
    dom % blocks(b) % weights      = 1.0
    dom % blocks(b) % face_weights = 1.0
  end do

  do b = 1, dom % n_blocks
    dom % blocks(b) % corners(0)=1       ! suppose it is properly oriented

    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(2),*) dom % blocks(b) % resolutions(1)
    read(line % tokens(3),*) dom % blocks(b) % resolutions(2)
    read(line % tokens(4),*) dom % blocks(b) % resolutions(3)

    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *)               &  ! block weights 
         dom % blocks(b) % weights(1),  &
         dom % blocks(b) % weights(2),  &
         dom % blocks(b) % weights(3)

    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *)                                             &
         dom % blocks(b) % corners(1), dom % blocks(b) % corners(2),  &
         dom % blocks(b) % corners(3), dom % blocks(b) % corners(4),  &
         dom % blocks(b) % corners(5), dom % blocks(b) % corners(6),  &
         dom % blocks(b) % corners(7), dom % blocks(b) % corners(8)

    !---------------------------!
    !   Check if the block is   ! 
    !     properly oriented     !
    !---------------------------!
    do n=1,8
      xt(n) = dom % points(dom % blocks(b) % corners(n)) % x
      yt(n) = dom % points(dom % blocks(b) % corners(n)) % y
      zt(n) = dom % points(dom % blocks(b) % corners(n)) % z
    end do

    if(Tet_Volume( xt(2),yt(2),zt(2), xt(5),yt(5),zt(5),  &
                   xt(3),yt(3),zt(3), xt(1),yt(1),zt(1) )  < 0) then
      dom % blocks(b) % corners(0)=-1            !  It's nor properly oriented
      call Swap_Integers(dom % blocks(b) % corners(2),  &
                         dom % blocks(b) % corners(3))
      call Swap_Integers(dom % blocks(b) % corners(6),  &
                         dom % blocks(b) % corners(7))
      call Swap_Reals(dom % blocks(b) % weights(1),  &
                      dom % blocks(b) % weights(2))
      dom % blocks(b) % weights(1) = 1.0 / dom % blocks(b) % weights(1)
      dom % blocks(b) % weights(2) = 1.0 / dom % blocks(b) % weights(2)
      call Swap_Integers(dom % blocks(b) % resolutions(1),  &
                         dom % blocks(b) % resolutions(2))
      print *, 'Warning: Block ',b,' was not properly oriented'
    end if
  end do                 ! through dom % blocks

  !-----------------------------!
  !   Set the corners of each   !
  !      face of the block      !
  !-----------------------------!
  do b = 1, dom % n_blocks
    do fc = 1,6
      do n = 1,4
        dom % blocks(b) % faces(fc, n) =  &
        dom % blocks(b) % corners(face_nodes(fc,n))
      end do
    end do
  end do

  !----------------------------------------------!
  !   Lines                                      !
  !----------------------------------------------!
  !   Lines can be prescribed point by point     ! 
  !   or with just a weighting factor.           !
  !----------------------------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) dom % n_lines  ! number of defined dom % lines

  call Domain_Mod_Allocate_Lines(dom, dom % n_lines)

  do l=1, dom % n_lines
    call Tokenizer_Mod_Read_Line(9)

    read(line % tokens(1),*) npnt
    read(line % tokens(2),*) dom % lines(l) % points(1)
    read(line % tokens(3),*) dom % lines(l) % points(2)

    call Find_Line(dom,                         &
                   dom % lines(l) % points(1),  &
                   dom % lines(l) % points(2),  &
                   dom % lines(l) % resolution)

    ! Does this need a more elegant solution?
    allocate(dom % lines(l) % x( dom % lines(l) % resolution ))
    allocate(dom % lines(l) % y( dom % lines(l) % resolution ))
    allocate(dom % lines(l) % z( dom % lines(l) % resolution ))

    ! Point by point
    if(npnt > 0) then
      do n=1,dom % lines(l) % resolution
        call Tokenizer_Mod_Read_Line(9)
        read(line % tokens(2),*) dom % lines(l) % x(n)
        read(line % tokens(3),*) dom % lines(l) % y(n)
        read(line % tokens(4),*) dom % lines(l) % z(n)
      end do

    ! Weight factor
    else
      call Tokenizer_Mod_Read_Line(9)
      read(line % tokens(1), *) dom % lines(l) % weight
    endif 

  end do

  !----------------------------------------!
  !   Copy block weights to face weights   !
  !----------------------------------------!
  do b = 1, dom % n_blocks
    do fc = 1,6                          !  face of the block
      dom % blocks(b) % face_weights(fc, 1) = dom % blocks(b) % weights(1)
      dom % blocks(b) % face_weights(fc, 2) = dom % blocks(b) % weights(2)
      dom % blocks(b) % face_weights(fc, 3) = dom % blocks(b) % weights(3)
    end do
  end do

  !--------------!
  !   Surfaces   !
  !--------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) nsurf     ! number of defined surfaces

  do s = 1, nsurf
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole,*) dum, n1, n2, n3, n4
    call Find_Surface(dom, n1, n2, n3, n4, b, fc)
    print *, '# block: ', b, ' surf: ', fc
    n = (b-1)*6 + fc         ! surface number

    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *) dom % blocks(b) % face_weights(fc,1),  &
                          dom % blocks(b) % face_weights(fc,2),  &
                          dom % blocks(b) % face_weights(fc,2)
  end do

  !---------------------------------------!
  !   Is there enough allocated memory?   !
  !---------------------------------------!

  ! Nodes & faces
  n_nodes_check = 0
  n_faces_check = 0
  do b = 1, dom % n_blocks
    ni = dom % blocks(b) % resolutions(1)
    nj = dom % blocks(b) % resolutions(2)
    nk = dom % blocks(b) % resolutions(3)
    n_nodes_check=n_nodes_check + ni*nj*nk
    n_faces_check=n_faces_check + ni*nj*nk + 2*( (ni*nj)+(nj*nk)+(ni*nk) )
  end do

  if( (n_faces_check > grid % max_n_faces) .or.  &
      (n_nodes_check > grid % max_n_nodes) ) then
    print *, '# Error message from TFlowS:'
  end if

  if( n_faces_check  > grid % max_n_faces ) then
    print *, '# The estimated number of sides is :', n_faces_check
    print *, '# There is space available only for:', grid % max_n_faces
    print *, '# Increase the parameter grid % max_n_faces in the input file'
    print *, '# and re-run the code !'
  end if

  if( n_nodes_check  > grid % max_n_nodes ) then
    print *, '# The estimated number of nodes is :', n_nodes_check
    print *, '# There is space available only for:', grid % max_n_nodes
    print *, '# Increase the parameter grid % max_n_nodes in the input file'
    print *, '# and re-run the code !'
  end if 

  if( (n_faces_check > grid % max_n_faces) .or.  &
      (n_nodes_check > grid % max_n_nodes) ) then
    stop
  end if

  !---------------------------------------!
  !   Boundary conditions and materials   !
  !---------------------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) dom % n_regions  ! number of regions (can be boundary 
                                             ! conditions or materials)

  call Domain_Mod_Allocate_Regions(dom, dom % n_regions)

  do n=1, dom % n_regions
    dom % regions(n) % face=''

    call Tokenizer_Mod_Read_Line(9)
    if(line % n_tokens .eq. 7) then
      read(line % whole,*)  dum,            &  
                   dom % regions(n) % is,   &
                   dom % regions(n) % js,   &
                   dom % regions(n) % ks,   & 
                   dom % regions(n) % ie,   &
                   dom % regions(n) % je,   &
                   dom % regions(n) % ke   
    else if(line % n_tokens .eq. 2) then
      read(line % tokens(1),*)       dum           
      read(line % tokens(2),'(A4)')  & 
           dom % regions(n) % face
      call To_Upper_Case(dom % regions(n) % face)           
    end if

    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(1), *) dom % regions(n) % block
    read(line % tokens(2), *) dom % regions(n) % name    
    call To_Upper_Case(dom % regions(n) % name)           

    ! if( dom % blocks(b_cond(n,7)) % points(0) .eq. -1 ) then
    !   call Swap_Integers( dom % regions(n) % is,dom % regions(n) % js )
    !   call Swap_Integers( dom % regions(n) % ie,dom % regions(n) % je )
    ! end if

  end do

  !-------------------------!
  !   Periodic boundaries   !
  !-------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *)  n_periodic_cond  ! number of periodic boundaries
  print '(a38,i7)', '# Number of periodic boundaries:     ', n_periodic_cond 

  allocate (periodic_cond(n_periodic_cond,8))

  do n=1,n_periodic_cond
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *) dum, periodic_cond(n,1), periodic_cond(n,2),  &
                               periodic_cond(n,3), periodic_cond(n,4) 
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *)      periodic_cond(n,5), periodic_cond(n,6),  &
                               periodic_cond(n,7), periodic_cond(n,8)
  end do

  !---------------------!
  !   Copy boundaries   !
  !---------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *)  n_copy_cond  ! number of copy boundaries
  print '(a38,i7)', '# Number of copy boundaries:         ', n_copy_cond

  allocate (copy_cond(n_copy_cond,8))

  do n=1,n_copy_cond
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *) dum, copy_cond(n,1), copy_cond(n,2),  &
                               copy_cond(n,3), copy_cond(n,4) 
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *)      copy_cond(n,5), copy_cond(n,6),  &
                               copy_cond(n,7),copy_cond(n,8)
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *)  copy_cond(n,0)
  end do

  !-----------------------------------!
  !   Refinement levels and regions   !
  !-----------------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) n_refine_levels ! number of refinement levels
  print '(a38,i7)', '# Number of refinement levels:       ', n_refine_levels

  allocate (refined_regions(n_refine_levels, 1024, 0:6))
  allocate (n_refined_regions(n_refine_levels))

  do l=1,n_refine_levels
    print *, 'Level: ', l
    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(2), *) n_refined_regions(l)  

    ! Browse through regions in level "l"
    do n = 1, n_refined_regions(l)
      call Tokenizer_Mod_Read_Line(9)
      read(line % tokens(3),*) answer
      call To_Upper_Case(answer)
      if(answer .eq. 'RECTANGLE') then
        refined_regions(l,n,0) = RECTANGLE
      elseif(answer .eq. 'ELIPSOID') then
        refined_regions(l,n,0) = ELIPSOID 
      elseif(answer .eq. 'PLANE') then
        refined_regions(l,n,0) = PLANE
      else
        print *, 'Error in input file: ', answer 
        stop
      endif 

      call Tokenizer_Mod_Read_Line(9)
      read(line % whole, *)                                       &
                refined_regions(l,n,1),refined_regions(l,n,2),    &
                refined_regions(l,n,3),refined_regions(l,n,4),    &
                refined_regions(l,n,5),refined_regions(l,n,6)   
    end do
  end do

  !-----------------------!
  !   Smoothing regions   !
  !-----------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) n_smoothing_regions  ! number of smoothing regions 

  print '(a38,i7)', '# Number of (non)smoothing regions:  ', n_smoothing_regions

  allocate (smooth_in_x   (n_smoothing_regions))
  allocate (smooth_in_y   (n_smoothing_regions))
  allocate (smooth_in_z   (n_smoothing_regions))
  allocate (smooth_iters  (n_smoothing_regions))
  allocate (smooth_relax  (n_smoothing_regions))
  allocate (smooth_regions(n_smoothing_regions, 0:6))

  do n=1, n_smoothing_regions
    smooth_in_x(n) = .false.
    smooth_in_y(n) = .false.
    smooth_in_z(n) = .false.
    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(1), *) smooth_regions(n,0)  
    if(line % n_tokens .eq. 4) then   ! smoothing in three directions
      smooth_in_x(n) = .true.
      smooth_in_y(n) = .true.
      smooth_in_z(n) = .true.
    else if(line % n_tokens .eq. 3) then
      call To_Upper_Case(line % tokens(2))
      call To_Upper_Case(line % tokens(3))
      if( line % tokens(2) .eq. 'X' )  &
          smooth_in_x(n) = .true.
      if( line % tokens(3) .eq. 'X' )  &
          smooth_in_x(n) = .true.
      if( line % tokens(2) .eq. 'Y' )  &
          smooth_in_y(n) = .true.
      if( line % tokens(3) .eq. 'Y' )  &
          smooth_in_y(n) = .true.
      if( line % tokens(2) .eq. 'Z' )  &
          smooth_in_z(n) = .true.
      if( line % tokens(3) .eq. 'Z' )  &
          smooth_in_z(n) = .true.
    else if(line % n_tokens .eq. 2) then
      call To_Upper_Case(line % tokens(2))
      if( line % tokens(2) .eq. 'X' )  &
          smooth_in_x(n) = .true.
      if( line % tokens(2) .eq. 'Y' )  &
          smooth_in_y(n) = .true.
      if( line % tokens(2) .eq. 'Z' )  &
          smooth_in_z(n) = .true.
    end if 

    ! Read the coordinates of the (non)smoothed region
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *) smooth_iters(n), smooth_relax(n)
    call Tokenizer_Mod_Read_Line(9)
    read(line % whole, *) smooth_regions(n,1),  &
                          smooth_regions(n,2),  &
                          smooth_regions(n,3),  &
                          smooth_regions(n,4),  &
                          smooth_regions(n,5),  &
                          smooth_regions(n,6)   
  end do

  close(9)

  end subroutine
