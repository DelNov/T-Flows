!==============================================================================!
  subroutine Sort_Cells_By_Thread(Grid)
!------------------------------------------------------------------------------!
!>  The subroutine is tailored for sorting cells in the computational grid
!>  based on the thread they are assigned to in an OpenMP environment. This
!>  sorting is crucial for efficient parallel computation
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Preparing for Sorting:                                                   !
!     - Allocates memory for temporary arrays and stores current cell data.    !
!     - Preserves old processor numbers for later validation.                  !
!   * Setting Sorting Criteria:                                                !
!     - Establishes the criteria for sorting based on threads and cell indices.!
!   * Performing the Sort:                                                     !
!     - Sorts cells according to the defined criteria, aligning them by thread.!
!   * Updating Grid Data:                                                      !
!     - Rearranges cell data, including topology and geometry, to reflect the  !
!       new sorted order.                                                      !
!   * Validating Post-Sort Structure:                                          !
!     - Ensures that the processor numbers remain consistent after sorting.    !
!                                                                              !
!   Note                                                                       !
!     The subroutine has a sister procedure in Sort_Cells_By_Coordinates,      !
!     which is a bad practice because there is a lot of code duplication, but  !
!     I will think about it at a later time.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! computational grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, c, c1, c2, n, nc
  integer              :: mcc, mcn, mcf
  integer, allocatable :: old_nn   (:)       ! old number of nodes (per cell)
  integer, allocatable :: old_nodes(:,:)     ! old nodes list per cell
  integer, allocatable :: old_nc   (:)       ! old number of cells (per cell)
  integer, allocatable :: old_nf   (:)       ! old number of faces (per cell)
  integer, allocatable :: old_faces(:,:)     ! old faces list per cell
  integer, allocatable :: old_cells(:,:)     ! old cells list per cell
  integer, allocatable :: thread(:), cell(:), old_proc(:)
!==============================================================================!

  Assert(PROGRAM_NAME .eq. 'Process')

  mcc = size(Grid % cells_c, 1)
  mcn = size(Grid % cells_n, 1)
  mcf = size(Grid % cells_f, 1)

  allocate(old_nn        (Grid % n_cells));  old_nn   (:)   = 0
  allocate(old_nodes(mcn, Grid % n_cells));  old_nodes(:,:) = 0
  allocate(old_nc        (Grid % n_cells));  old_nc   (:)   = 0
  allocate(old_cells(mcc ,Grid % n_cells));  old_cells(:,:) = 0
  allocate(old_nf        (Grid % n_cells));  old_nf   (:)   = 0
  allocate(old_faces(mcf, Grid % n_cells));  old_faces(:,:) = 0
  allocate(thread        (Grid % n_cells));  thread   (:)   = 0
  allocate(cell          (Grid % n_cells));  cell     (:)   = 0
  allocate(old_proc      (Grid % n_cells))

  !---------------------------------!
  !   Store old processor numbers   !
  !---------------------------------!
  do c = 1, Grid % n_cells
    old_proc(c) = Grid % Comm % cell_proc(c)
  end do

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, Grid % n_cells
    old_nn   (c)        = Grid % cells_n_nodes(c)
    old_nodes(1:mcn, c) = Grid % cells_n(1:mcn, c)
    old_nf   (c)        = Grid % cells_n_faces(c)
    old_faces(1:mcf, c) = Grid % cells_f(1:mcf, c)
    old_nc   (c)        = Grid % cells_n_cells(c)
    old_cells(1:mcc, c) = Grid % cells_c(1:mcc, c)
  end do

  !--------------------------!
  !   Set sorting criteria   !
  !--------------------------!
  do c = 1, Grid % n_cells
    thread(c) = Grid % Omp % cell_thread(c)  ! in buffers, must be > n_threads
    cell  (c) = c
    Grid % old_c(c) = c
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Two_Int_Carry_Int(thread(1:Grid % n_cells),  &
                                cell  (1:Grid % n_cells),  &
                                Grid % old_c(1:Grid % n_cells))

  ! Check if buffers stayed untouched
  do c = Cells_In_Buffers()
    Assert(Grid % old_c(c) .eq. c)
  end do

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
    Grid % cells_n_nodes(Grid % new_c(c)) = old_nn          (c)
    Grid % cells_n(1:mcn,Grid % new_c(c)) = old_nodes(1:mcn, c)
    Grid % cells_n_faces(Grid % new_c(c)) = old_nf          (c)
    Grid % cells_f(1:mcf,Grid % new_c(c)) = old_faces(1:mcf, c)
    Grid % cells_n_cells(Grid % new_c(c)) = old_nc          (c)
    Grid % cells_c(1:mcc,Grid % new_c(c)) = old_cells(1:mcc, c)
  end do
  nc = Grid % n_cells  ! abbreviate the syntax
  call Sort % Real_By_Index(nc, Grid % xc                (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % yc                (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % zc                (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % vol               (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % ixx               (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % iyy               (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % izz               (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % ixy               (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % ixz               (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % iyz               (1), Grid % new_c(1))
  call Sort % Real_By_Index(nc, Grid % wall_dist         (1), Grid % new_c(1))
  call Sort % Int_By_Index (nc, Grid % por               (1), Grid % new_c(1))
  ! Sorting of the following three fields does not exist in the sister function
  call Sort % Int_By_Index (nc, Grid % Omp  % cell_thread(1), Grid % new_c(1))
  call Sort % Int_By_Index (nc, Grid % Comm % cell_glo   (1), Grid % new_c(1))
  call Sort % Int_By_Index (nc, Grid % Comm % cell_proc  (1), Grid % new_c(1))

  !-----------------------------!
  !   Check processor numbers   !
  !-----------------------------!
  do c = 1, Grid % n_cells
    Assert(old_proc(c) .eq. Grid % Comm % cell_proc(c))
  end do

  end subroutine
