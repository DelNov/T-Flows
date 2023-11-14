!==============================================================================!
  subroutine Load_Dom(Generate, Dom, smr, ref, Grid)
!------------------------------------------------------------------------------!
!>  Routine for loading and initializing various aspects of a computational
!>  domain from a domain file (.dom file).
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!   * Initialization and Memory Allocation: It starts by allocating memory for !
!     different grid components such as nodes, cells, and faces. This includes !
!     setting up the grid size and other basic parameters.                     !
!   * Reading Domain File: The core functionality revolves around reading the  !
!     .dom file, which contains detailed specifications of the grid. This      !
!     includes the geometry of the domain (corners, blocks, lines), boundary   !
!     conditions, refinement levels, and smoothing ranges.                     !
!   * Setting Up Geometrical Elements:                                         !
!     - Corners and Blocks: It reads and processes the coordinates of domain   !
!       points (corners) and the specifications of grid blocks.                !
!     - Lines: The subroutine handles the definition of lines, which can be    !
!       specified point by point or with weighting factors.                    !
!     - Surfaces: It processes the surface elements of the blocks, ensuring    !
!       proper orientation and weighting.                                      !
!   * Boundary Conditions and Materials: It reads ranges for boundary          !
!     conditions and materials.                                                !
!   * Periodic Boundaries: The subroutine handles periodic boundary conditions !
!   * Refinement and Smoothing Ranges: It sets up refinement levels and        !
!     smoothing ranges as specified in the domain file.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Generate_Type) :: Generate  !! parent class
  type(Domain_Type)    :: Dom       !! computational domain
  type(Smooths_Type)   :: smr       !! smoothing regions
  type(Refines_Type)   :: ref       !! refinement regions
  type(Grid_Type)      :: Grid      !! computational grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: b, i, l, s, i_fac, n, n1, n2, n3, n4, dumi, fu
  integer       :: n_faces_check, n_nodes_check
  integer       :: ni, nj, nk, npnt, nsurf
  character(SL) :: dum
  character(SL) :: domain_name
  character(SL) :: answer
  real          :: xt(8), yt(8), zt(8)
  integer       :: fn(6,4)
!------------------------------------------------------------------------------!
  include 'Block_Numbering.h90'
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Generate)
!==============================================================================!

  ! Copy face-node numbering for blocks
  fn = hex_block

  !---------------------------------------------!
  !   The Grid will not have polyhedral cells   !
  !---------------------------------------------!
  Grid % polyhedral = .false.

  print '(a)', ' #========================================'
  print '(a)', ' # Input problem name: (without extension)'
  print '(a)', ' #----------------------------------------'
  call File % Read_Line(5)
  read(Line % tokens(1), *) problem_name(1)

  Grid % name = problem_name(1)
  call String % To_Upper_Case(Grid % name)

  call File % Set_Name(domain_name, extension='.dom')
  call File % Open_For_Reading_Ascii(domain_name, fu)

  !-----------------------------------------------------------------!
  !   Max. number of nodes (cells), boundary faces and cell faces   !
  !-----------------------------------------------------------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) Grid % max_n_nodes
  read(Line % tokens(2), *) Grid % max_n_bnd_cells
  read(Line % tokens(3), *) Grid % max_n_faces

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  print '(a)',      ' # Allocating memory for: '
  print '(a,i9,a)', ' #', Grid % max_n_nodes,     ' nodes and cells'
  print '(a,i9,a)', ' #', Grid % max_n_bnd_cells, ' boundary cells'
  print '(a,i9,a)', ' #', Grid % max_n_faces,     ' cell faces'

  allocate (Grid % region % at_cell(-Grid % max_n_bnd_cells-1:-1))
  Grid % region % at_cell = 0

  ! Variables in Grid_Mod
  call Grid % Allocate_Nodes(Grid % max_n_nodes)

  call Grid % Allocate_Cells(Grid % max_n_nodes,      &
                             Grid % max_n_bnd_cells)

  call Grid % Allocate_Faces(Grid % max_n_faces, 0)

  call Refines_Mod_Allocate_Cells(ref,                     &
                                  Grid % max_n_bnd_cells,  &
                                  Grid % max_n_nodes)

  ! Variables still declared in Gen_Mod.h90:
  allocate (face_c_to_c(Grid % max_n_faces,2))
  face_c_to_c = 0
  allocate (twin_n(Grid % max_n_nodes,0:8))
  twin_n (:,:) = 0

  print *, '# Allocation successfull !'

  !-------------!
  !   Corners   !
  !-------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) Dom % n_points  ! number of points

  call Dom % Allocate_Points(Dom % n_points)

  do i = 1, Dom % n_points
    call File % Read_Line(fu)
    read(Line % tokens(2),*) Dom % points(i) % x
    read(Line % tokens(3),*) Dom % points(i) % y
    read(Line % tokens(4),*) Dom % points(i) % z
  end do

  !------------!
  !   Blocks   !
  !------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) Dom % n_blocks  ! number of blocks

  call Dom % Allocate_Blocks(Dom % n_blocks)

  ! Initialize weights
  do b=1, Dom % n_blocks
    Dom % blocks(b) % weights      = 1.0
    Dom % blocks(b) % face_weights = 1.0
  end do

  do b = 1, Dom % n_blocks
    Dom % blocks(b) % corners(0)=1       ! suppose it is properly oriented

    call File % Read_Line(fu)
    read(Line % tokens(2),*) Dom % blocks(b) % resolutions(1)
    read(Line % tokens(3),*) Dom % blocks(b) % resolutions(2)
    read(Line % tokens(4),*) Dom % blocks(b) % resolutions(3)

    call File % Read_Line(fu)
    read(Line % whole, *)               &  ! block weights
         Dom % blocks(b) % weights(1),  &
         Dom % blocks(b) % weights(2),  &
         Dom % blocks(b) % weights(3)

    call File % Read_Line(fu)
    read(Line % whole, *)                                             &
         Dom % blocks(b) % corners(1), Dom % blocks(b) % corners(2),  &
         Dom % blocks(b) % corners(3), Dom % blocks(b) % corners(4),  &
         Dom % blocks(b) % corners(5), Dom % blocks(b) % corners(6),  &
         Dom % blocks(b) % corners(7), Dom % blocks(b) % corners(8)

    !---------------------------!
    !   Check if the block is   !
    !     properly oriented     !
    !---------------------------!
    do n=1,8
      xt(n) = Dom % points(Dom % blocks(b) % corners(n)) % x
      yt(n) = Dom % points(Dom % blocks(b) % corners(n)) % y
      zt(n) = Dom % points(Dom % blocks(b) % corners(n)) % z
    end do

    if(Math % Tet_Volume( xt(2),yt(2),zt(2), xt(5),yt(5),zt(5),  &
                          xt(3),yt(3),zt(3), xt(1),yt(1),zt(1) )  < 0) then
      Dom % blocks(b) % corners(0)=-1            !  It's nor properly oriented
      call Swap_Int(Dom % blocks(b) % corners(2),  &
                    Dom % blocks(b) % corners(3))
      call Swap_Int(Dom % blocks(b) % corners(6),  &
                    Dom % blocks(b) % corners(7))
      call Swap_Real(Dom % blocks(b) % weights(1),  &
                     Dom % blocks(b) % weights(2))
      Dom % blocks(b) % weights(1) = 1.0 / Dom % blocks(b) % weights(1)
      Dom % blocks(b) % weights(2) = 1.0 / Dom % blocks(b) % weights(2)
      call Swap_Int(Dom % blocks(b) % resolutions(1),  &
                    Dom % blocks(b) % resolutions(2))
      print *, 'Warning: Block ',b,' was not properly oriented'
    end if
  end do                 ! through Dom % blocks

  !-----------------------------!
  !   Set the corners of each   !
  !      face of the block      !
  !-----------------------------!
  do b = 1, Dom % n_blocks
    do i_fac = 1, 6
      do n = 1, 4
        Dom % blocks(b) % faces(i_fac, n) =  &
        Dom % blocks(b) % corners(fn(i_fac,n))
      end do
    end do
  end do

  !----------------------------------------------!
  !   Lines                                      !
  !----------------------------------------------!
  !   Lines can be prescribed point by point     !
  !   or with just a weighting factor.           !
  !----------------------------------------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) Dom % n_lines  ! number of defined Dom % lines

  call Dom % Allocate_Lines(Dom % n_lines)

  do l=1, Dom % n_lines
    call File % Read_Line(fu)

    read(Line % tokens(1),*) npnt
    read(Line % tokens(2),*) Dom % lines(l) % points(1)
    read(Line % tokens(3),*) Dom % lines(l) % points(2)

    call Dom % Find_Line(Dom % lines(l) % points(1),  &
                         Dom % lines(l) % points(2),  &
                         Dom % lines(l) % resolution)

    ! Does this need a more elegant solution?
    allocate(Dom % lines(l) % x( Dom % lines(l) % resolution ))
    allocate(Dom % lines(l) % y( Dom % lines(l) % resolution ))
    allocate(Dom % lines(l) % z( Dom % lines(l) % resolution ))

    ! Point by point
    if(npnt > 0) then
      do n=1,Dom % lines(l) % resolution
        call File % Read_Line(fu)
        read(Line % tokens(2),*) Dom % lines(l) % x(n)
        read(Line % tokens(3),*) Dom % lines(l) % y(n)
        read(Line % tokens(4),*) Dom % lines(l) % z(n)
      end do

    ! Weight factor
    else
      call File % Read_Line(fu)
      read(Line % tokens(1), *) Dom % lines(l) % weight
    end if

  end do

  !----------------------------------------!
  !   Copy block weights to face weights   !
  !----------------------------------------!
  do b = 1, Dom % n_blocks
    do i_fac = 1,6                          !  face of the block
      Dom % blocks(b) % face_weights(i_fac, 1) = Dom % blocks(b) % weights(1)
      Dom % blocks(b) % face_weights(i_fac, 2) = Dom % blocks(b) % weights(2)
      Dom % blocks(b) % face_weights(i_fac, 3) = Dom % blocks(b) % weights(3)
    end do
  end do

  !--------------!
  !   Surfaces   !
  !--------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) nsurf     ! number of defined surfaces

  do s = 1, nsurf
    call File % Read_Line(fu)
    read(Line % whole,*) dum, n1, n2, n3, n4
    call Dom % Find_Surface(n1, n2, n3, n4, b, i_fac)
    print *, '# block: ', b, ' surf: ', i_fac
    n = (b-1)*6 + i_fac         ! surface number

    call File % Read_Line(fu)
    read(Line % whole, *) Dom % blocks(b) % face_weights(i_fac,1),  &
                          Dom % blocks(b) % face_weights(i_fac,2),  &
                          Dom % blocks(b) % face_weights(i_fac,2)
  end do

  !---------------------------------------!
  !   Is there enough allocated memory?   !
  !---------------------------------------!

  ! Nodes & faces
  n_nodes_check = 0
  n_faces_check = 0
  do b = 1, Dom % n_blocks
    ni = Dom % blocks(b) % resolutions(1)
    nj = Dom % blocks(b) % resolutions(2)
    nk = Dom % blocks(b) % resolutions(3)
    n_nodes_check=n_nodes_check + ni*nj*nk
    n_faces_check=n_faces_check + ni*nj*nk + 2*( (ni*nj)+(nj*nk)+(ni*nk) )
  end do

  if( (n_faces_check > Grid % max_n_faces) .or.  &
      (n_nodes_check > Grid % max_n_nodes) ) then
    print *, '# Error message from T-Flows:'
  end if

  if( n_faces_check  > Grid % max_n_faces ) then
    print *, '# The estimated number of faces is :', n_faces_check
    print *, '# There is space available only for:', Grid % max_n_faces
    print *, '# Increase the parameter Grid % max_n_faces in the input file'
    print *, '# and re-run the code !'
  end if

  if( n_nodes_check  > Grid % max_n_nodes ) then
    print *, '# The estimated number of nodes is :', n_nodes_check
    print *, '# There is space available only for:', Grid % max_n_nodes
    print *, '# Increase the parameter Grid % max_n_nodes in the input file'
    print *, '# and re-run the code !'
  end if

  if( (n_faces_check > Grid % max_n_faces) .or.  &
      (n_nodes_check > Grid % max_n_nodes) ) then
    stop
  end if

  !---------------------------------------!
  !   Boundary conditions and materials   !
  !---------------------------------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) Dom % n_ranges   ! number of ranges (can be bnd.
                                             ! conditions or materials)
  call Dom % Allocate_Ranges(Dom % n_ranges)

  do n = 1, Dom % n_ranges
    Dom % ranges(n) % face=''

    call File % Read_Line(fu)
    if(Line % n_tokens .eq. 7) then
      read(Line % whole,*)  dum,           &
                   Dom % ranges(n) % is,   &
                   Dom % ranges(n) % js,   &
                   Dom % ranges(n) % ks,   &
                   Dom % ranges(n) % ie,   &
                   Dom % ranges(n) % je,   &
                   Dom % ranges(n) % ke
    else if(Line % n_tokens .eq. 2) then
      read(Line % tokens(1),*)       dum
      read(Line % tokens(2),'(A4)')  &
           Dom % ranges(n) % face
      call String % To_Upper_Case(Dom % ranges(n) % face)
    end if

    call File % Read_Line(fu)
    read(Line % tokens(1), *) Dom % ranges(n) % block
    read(Line % tokens(2), *) Dom % ranges(n) % name
    call String % To_Upper_Case(Dom % ranges(n) % name)

    ! if( Dom % blocks(b_cond(n,7)) % points(0) .eq. -1 ) then
    !   call Swap_Int( Dom % ranges(n) % is,Dom % ranges(n) % js )
    !   call Swap_Int( Dom % ranges(n) % ie,Dom % ranges(n) % je )
    ! end if

  end do

  !-------------------------!
  !   Periodic boundaries   !
  !-------------------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *)  n_periodic_cond  ! number of periodic boundaries
  print '(a38,i7)', '# Number of periodic boundaries:     ', n_periodic_cond

  allocate (periodic_cond(n_periodic_cond,8))

  do n=1,n_periodic_cond
    call File % Read_Line(fu)
    read(Line % whole, *) dum, periodic_cond(n,1), periodic_cond(n,2),  &
                               periodic_cond(n,3), periodic_cond(n,4)
    call File % Read_Line(fu)
    read(Line % whole, *)      periodic_cond(n,5), periodic_cond(n,6),  &
                               periodic_cond(n,7), periodic_cond(n,8)
  end do

  !---------------------!
  !   Copy boundaries   !
  !---------------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *)  dumi  ! used to be number of copy boundaries

  !----------------------------------!
  !   Refinement levels and ranges   !
  !----------------------------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) ref % n_levels     ! number of refinement levels
  print '(a38,i7)', '# Number of refinement levels:       ', ref % n_levels

  call Refines_Mod_Allocate_Levels(ref, ref % n_levels)

  do l = 1, ref % n_levels
    print *, '# Level: ', l
    call File % Read_Line(fu)
    read(Line % tokens(2), *) ref % n_ranges(l)

    ! Browse through ranges in level "l"
    do n = 1, ref % n_ranges(l)
      call File % Read_Line(fu)
      read(Line % tokens(3),*) answer
      call String % To_Upper_Case(answer)
      ref % range(l,n) % shape = -1
      if(answer .eq. 'RECTANGLE') ref % range(l,n) % shape = RECTANGLE
      if(answer .eq. 'ELIPSOID')  ref % range(l,n) % shape = ELIPSOID
      if(answer .eq. 'PLANE')     ref % range(l,n) % shape = PLANE
      if(ref % range(l,n) % shape .eq. -1) then
        print *, 'ERROR!  Refinement shape not specified well by: ', answer
        stop
      end if

      call File % Read_Line(fu)
      read(Line % whole, *)                     &
                ref % range(l,n) % pnt(1) % x,  &
                ref % range(l,n) % pnt(1) % y,  &
                ref % range(l,n) % pnt(1) % z,  &
                ref % range(l,n) % pnt(2) % x,  &
                ref % range(l,n) % pnt(2) % y,  &
                ref % range(l,n) % pnt(2) % z
    end do
  end do

  !----------------------!
  !   Smoothing ranges   !
  !----------------------!
  call File % Read_Line(fu)
  read(Line % tokens(1), *) smr % n_smooths  ! number of smoothing ranges

  print '(a38,i7)', '# Number of (non)smoothing ranges:   ', smr % n_smooths

  call Smooths_Mod_Allocate_Smooths(smr, smr % n_smooths)

  do n = 1, smr % n_smooths
    smr % in_x(n) = .false.
    smr % in_y(n) = .false.
    smr % in_z(n) = .false.
    call File % Read_Line(fu)
    read(Line % tokens(1), *) dumi    ! this Line is probably not needed
    if(Line % n_tokens .eq. 4) then   ! smoothing in three directions
      smr % in_x(n) = .true.
      smr % in_y(n) = .true.
      smr % in_z(n) = .true.
    else if(Line % n_tokens .eq. 3) then
      call String % To_Upper_Case(Line % tokens(2))
      call String % To_Upper_Case(Line % tokens(3))
      if( Line % tokens(2) .eq. 'X' )  smr % in_x(n) = .true.
      if( Line % tokens(3) .eq. 'X' )  smr % in_x(n) = .true.
      if( Line % tokens(2) .eq. 'Y' )  smr % in_y(n) = .true.
      if( Line % tokens(3) .eq. 'Y' )  smr % in_y(n) = .true.
      if( Line % tokens(2) .eq. 'Z' )  smr % in_z(n) = .true.
      if( Line % tokens(3) .eq. 'Z' )  smr % in_z(n) = .true.
    else if(Line % n_tokens .eq. 2) then
      call String % To_Upper_Case(Line % tokens(2))
      if( Line % tokens(2) .eq. 'X' )  smr % in_x(n) = .true.
      if( Line % tokens(2) .eq. 'Y' )  smr % in_y(n) = .true.
      if( Line % tokens(2) .eq. 'Z' )  smr % in_z(n) = .true.
    end if

    ! Read the coordinates of the (non)smoothed region
    call File % Read_Line(fu)
    read(Line % whole, *) smr % iters(n),  &
                          smr % relax(n)

    call File % Read_Line(fu)
    read(Line % whole, *) smr % x_min(n),   &
                          smr % y_min(n),   &
                          smr % z_min(n),   &
                          smr % x_max(n),   &
                          smr % y_max(n),   &
                          smr % z_max(n)
  end do

  close(fu)

  end subroutine
