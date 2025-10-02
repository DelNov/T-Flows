!==============================================================================!
  subroutine Sort_Cells_By_Coordinates(Grid)
!------------------------------------------------------------------------------!
!>  Reorders cells based on their geometric positions. This subroutine is
!>  called from main functions within the Generate and Convert sub-programs,
!>  just after the calculation of geometric quantities of the cells.
!>  The subroutine hypothesizes that, since many computational domain in CFD
!>  are elongated in streamwise direction, sorting cells in that same
!>  direction, might lead to reduced bandwidth of matrices.  Although I
!>  didn't see any clear evidence of that in practical simulations, it is
!>  probably good to keep it around because its functionality might be put
!>  in use at a certain point in the future.
!------------------------------------------------------------------------------!
!   * Storing original cell data: The subroutine stores the original data of   !
!     the cells, including the number of nodes per cell, node indices, face    !
!     indices, neighboring cells, and boundary regions. This is crucial for    !
!     maintaining the integrity of the cell data during the sorting process.   !
!   * Sorting operation: The subroutine sorts the cells based on their center  !
!     coordinates (x, y, z). This is done with Sort % Three_Real_Carry_Int     !
!     function, which orders the cells based on their 3D center coordinates    !
!     while carrying along their original indices for reference.               !
!   * Updating references: The subroutine updates references in various grid   !
!     components, including faces and neighboring cells, to correspond to the  !
!     new cell ordering. This involves adjusting indices in face-to-cell       !
!     connections, neighboring cell info, and other geometric quantities.      !
!   * Calculating maximum differences: After reordering, the subroutine        !
!     calculates the maximum difference in cell indices between adjacent       !
!     faces. This provides an indication of the extent of reordering that      !
!     occurred, which can be important for understanding the reordering's      !
!     impact on the grid's structure.                                          !
!   * Conclusion and output: Finally, the subroutine prints information about  !
!     the maximum differences in cell indices, providing a summary of the      !
!     sorting's effect on the grid.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, c, c1, c2, n, nc
  integer              :: mcc, mcn, mcf, mcb
  integer              :: max_diff_1, max_diff_2, c1_s1, c2_s1, c1_s2, c2_s2
  integer, allocatable :: old_nn   (:)       ! old number of nodes (per cell)
  integer, allocatable :: old_nodes(:,:)     ! old nodes list per cell
  integer, allocatable :: old_nc   (:)       ! old number of cells (per cell)
  integer, allocatable :: old_nf   (:)       ! old number of faces (per cell)
  integer, allocatable :: old_faces(:,:)     ! old faces list per cell
  integer, allocatable :: old_cells(:,:)     ! old cells list per cell
  integer, allocatable :: old_b_reg(:,:)     ! old cells boundary regions
!==============================================================================!

  Assert(PROGRAM_NAME .eq. 'Generate' .or. PROGRAM_NAME .eq. 'Convert')

  if(Math % Approx_Real(maxval(Grid % xc(:)), 0.0) .and.  &
     Math % Approx_Real(maxval(Grid % yc(:)), 0.0) .and.  &
     Math % Approx_Real(maxval(Grid % zc(:)), 0.0) .and.  &
     Math % Approx_Real(minval(Grid % xc(:)), 0.0) .and.  &
     Math % Approx_Real(minval(Grid % yc(:)), 0.0) .and.  &
     Math % Approx_Real(minval(Grid % zc(:)), 0.0)) then
    call Message % Error(80,                                                  &
            'You called the function Grid % Sort_Cells_By_Coordinates   ' //  &
            'before the cell coordinates have been calculated.  The     ' //  &
            'code can''t continue like this. \n Something is wrong with ' //  &
            'the logic of the algorithm, you better re-assess it.',           &
            file=__FILE__, line=__LINE__)
  end if

  mcc = size(Grid % cells_c,          1)
  mcn = size(Grid % cells_n,          1)
  mcf = size(Grid % cells_f,          1)
  mcb = size(Grid % cells_bnd_region, 1)

  allocate(old_nn        (Grid % n_cells));  old_nn   (:)   = 0
  allocate(old_nodes(mcn, Grid % n_cells));  old_nodes(:,:) = 0
  allocate(old_nc        (Grid % n_cells));  old_nc   (:)   = 0
  allocate(old_cells(mcc ,Grid % n_cells));  old_cells(:,:) = 0
  allocate(old_nf        (Grid % n_cells));  old_nf   (:)   = 0
  allocate(old_faces(mcf, Grid % n_cells));  old_faces(:,:) = 0
  allocate(old_b_reg(mcb, Grid % n_cells));  old_b_reg(:,:) = 0

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, Grid % n_cells
    old_nn   (c)        = Grid % cells_n_nodes          (c)
    old_nodes(1:mcn, c) = Grid % cells_n         (1:mcn, c)
    old_nf   (c)        = Grid % cells_n_faces          (c)
    old_faces(1:mcf, c) = Grid % cells_f         (1:mcf, c)
    old_nc   (c)        = Grid % cells_n_cells          (c)
    old_cells(1:mcc, c) = Grid % cells_c         (1:mcc, c)
    old_b_reg(1:mcb, c) = Grid % cells_bnd_region(1:mcb, c)
  end do

  !----------------------------!
  !   Store old cell numbers   !
  !----------------------------!
  do c = 1, Grid % n_cells
    Grid % old_c(c) = c
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Three_Real_Carry_Int(Grid % xc(1:Grid % n_cells),  &
                                   Grid % yc(1:Grid % n_cells),  &
                                   Grid % zc(1:Grid % n_cells),  &
                                   Grid % old_c(1:Grid % n_cells))
  ! This is a bit of a bluff
  do c = 1, Grid % n_cells
    Grid % new_c(Grid % old_c(c)) = c
  end do

  !----------------------------------!
  !   Update cell numbers at faces   !
  !----------------------------------!
  do s = 1, Grid % n_faces + Grid % n_shadows
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    Grid % faces_c(1, s) = Grid % new_c(c1)
    if(c2 > 0) then
      Grid % faces_c(2, s) = Grid % new_c(c2)

      ! If the face changed its orientation during cell renumeration
      if(Grid % faces_c(2, s) < Grid % faces_c(1, s)) then

        ! Swap c1 and c2 of course ...
        call Swap_Int(Grid % faces_c(1, s), Grid % faces_c(2, s))

        ! ...but also reverse the order of face's nodes ...
        n = Grid % faces_n_nodes(s)  ! number of nodes in this face
        call Sort % Reverse_Order_Int(Grid % faces_n(1:n, s))

        ! ... and fix the geometrical quantities
        Grid % sx(s) = -Grid % sx(s)
        Grid % sy(s) = -Grid % sy(s)
        Grid % sz(s) = -Grid % sz(s)
        Grid % dx(s) = -Grid % dx(s)
        Grid % dy(s) = -Grid % dy(s)
        Grid % dz(s) = -Grid % dz(s)
        Grid % f (s) = 1.0 - Grid % f(s)
      end if
    end if
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = 1, Grid % n_cells
    Grid % cells_n_nodes         (Grid % new_c(c)) = old_nn          (c)
    Grid % cells_n         (1:mcn,Grid % new_c(c)) = old_nodes(1:mcn, c)
    Grid % cells_n_faces         (Grid % new_c(c)) = old_nf          (c)
    Grid % cells_f         (1:mcf,Grid % new_c(c)) = old_faces(1:mcf, c)
    Grid % cells_n_cells         (Grid % new_c(c)) = old_nc          (c)
    Grid % cells_c         (1:mcc,Grid % new_c(c)) = old_cells(1:mcc, c)
    Grid % cells_bnd_region(1:mcb,Grid % new_c(c)) = old_b_reg(1:mcb, c)
  end do
  nc = Grid % n_cells  ! abbreviate the syntax
  call Sort % Real_By_Index(nc, Grid % vol      (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % ixx      (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % iyy      (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % izz      (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % ixy      (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % ixz      (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % iyz      (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % wall_dist(1), Grid % new_c(1))
  call Sort % Int_By_Index (nc, Grid % por      (1), Grid % new_c(1))

  !--------------------------------------------!
  !   Find out distance between cell indices   !
  !--------------------------------------------!
  max_diff_1 = 0
  max_diff_2 = 0
  do s = 1, Grid % n_faces  ! + Grid % n_shadows
    c1_s1 = Grid % faces_c(1, s)
    c2_s1 = Grid % faces_c(2, s)
    if(c2_s1 > 0) then
      max_diff_1 = max((c2_s1 - c1_s1), max_diff_1)
      if(s < Grid % n_faces) then
        c1_s2 = Grid % faces_c(1, s+1)
        c2_s2 = Grid % faces_c(2, s+1)
        max_diff_2 = max(abs(c1_s1 - c1_s2), max_diff_2)
        if(c2_s1 > 0) then
          max_diff_2 = max(abs(c2_s1 - c2_s2), max_diff_2)
        end if
      end if
    end if
  end do
  print '(a)',    ' #=========================================================='
  print '(a)',    ' # In Sort_Cells_By_Coordinates'
  print '(a)',    ' #----------------------------------------------------------'
  print '(a,i9)', ' # Maximum cell difference at single face:   ', max_diff_1
  print '(a,i9)', ' # Maximum cell difference betwen two faces: ', max_diff_2

  end subroutine
