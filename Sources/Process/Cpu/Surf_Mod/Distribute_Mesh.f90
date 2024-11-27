!==============================================================================!
  subroutine Distribute_Mesh(Surf, verbose)
!------------------------------------------------------------------------------!
!>  This subroutine is essential for distributing the surface data across all
!>  processors in parallel simulations. Unlike Front_Type, which is divided
!>  among processors, Surf_Type shares the entire surface mesh among all
!>  processors. This subroutine synchronizes the mesh data, ensuring that each
!>  processor has the same mesh information, facilitating consistent
!>  computation across the parallel environment.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization and cleaning buffers: Sets up aliases for easier          !
!     reference and cleans buffers to prepare for data distribution.           !
!   * Estimation of mesh size per processor: Calculates the number of          !
!     vertices and elements each processor will handle and aggregates this     !
!     information globally.                                                    !
!   * Verbose output: Optionally prints cumulative numbers of elements and     !
!     vertices found, useful for debugging and verification in verbose mode.   !
!   * Accumulation calculation: Determines the cumulative number of vertices   !
!     and elements each processor is responsible for up to the current.        !
!   * Exchange of vertex coordinates: Populates buffers with vertex            !
!     coordinates and exchanges this information among processors.             !
!   * Summation of buffer arrays: Aggregates the buffered vertex coordinates   !
!     and element vertices across all processors.                              !
!   * Final assignment of coordinates and elements: Assigns the aggregated     !
!     vertex coordinates and elements back to the respective arrays, ensuring  !
!     each processor has complete and updated mesh information.                !
!   * Mesh size update: Updates the total number of vertices and elements in   !
!     the Surf_Type to reflect the complete mesh.                              !
!   * Ensuring consistency: After execution, all processors have complete and  !
!     consistent information about the entire mesh, crucial for the integrity  !
!     of parallel simulations.                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type),  target :: Surf     !! parent class
  logical                   :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Vert_Type),   pointer :: Vert(:)
  type(Elem_Type),   pointer :: Elem(:)
  integer,           pointer :: nv, ne
  integer                    :: j
  integer                    :: nv_tot, ne_tot, n_acc_vert, n_acc_elem
  integer, allocatable       :: nv_proc(:), ne_proc(:)
!==============================================================================!

  ! Take aliases
  Grid => Surf % pnt_grid
  Flow => Surf % pnt_flow
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  Elem => Surf % Elem

  ! Clean buffers
  Surf % buff_x(:) = 0.0
  Surf % buff_y(:) = 0.0
  Surf % buff_z(:) = 0.0
  Surf % buff_v(:) = 0.0
  Surf % buff_i(:) = 0
  Surf % buff_j(:) = 0
  Surf % buff_k(:) = 0
  Surf % buff_n(:) = 0

  !------------------------------------------------------------!
  !   Estimate number of vertices and elements per processor   !
  !------------------------------------------------------------!
  allocate(nv_proc(N_Procs())); nv_proc(:) = 0
  allocate(ne_proc(N_Procs())); ne_proc(:) = 0

  nv_proc(This_Proc()) = nv
  ne_proc(This_Proc()) = ne

  call Global % Sum_Int_Array(N_Procs(), nv_proc)
  call Global % Sum_Int_Array(N_Procs(), ne_proc)

  nv_tot = sum(nv_proc(1:N_Procs()))
  ne_tot = sum(ne_proc(1:N_Procs()))

  if(verbose) then
    if(First_Proc()) then
      print '(a40,i8)', ' # Cummulative number of elements found:', ne_tot
      print '(a40,i8)', ' # Cummulative number of vertices found:', nv_tot
    end if
  end if

  !-------------------------------------------------------!
  !   Estimate accumulated number of nodes and elements   !
  !-------------------------------------------------------!
  n_acc_vert = 0  ! accumulated nodes
  n_acc_elem = 0  ! accumulated nodes

  n_acc_vert = sum(nv_proc(1:This_Proc()-1))
  n_acc_elem = sum(ne_proc(1:This_Proc()-1))

  ! PRINT *, 'N_ACC_VERT', THIS_PROC, N_ACC_VERT
  ! PRINT *, 'N_ACC_ELEM', THIS_PROC, N_ACC_ELEM

  !------------------------------------------------------!
  !   Exchange vertex coordinates among the processors   !
  !------------------------------------------------------!

  do j = 1, nv
    Surf % buff_x(n_acc_vert+j) = Vert(j) % x_n
    Surf % buff_y(n_acc_vert+j) = Vert(j) % y_n
    Surf % buff_z(n_acc_vert+j) = Vert(j) % z_n
    ! WRITE(100+THIS_PROC, *)  N_ACC_VERT+J,  &
    !                          VERT(J) % X_N, VERT(J) % Y_N, VERT(J) % Z_N
  end do

  do j = 1, ne
    Surf % buff_i(n_acc_elem+j) = Elem(j) % v(1) + n_acc_vert
    Surf % buff_j(n_acc_elem+j) = Elem(j) % v(2) + n_acc_vert
    Surf % buff_k(n_acc_elem+j) = Elem(j) % v(3) + n_acc_vert
    Surf % buff_n(n_acc_elem+j) = Elem(j) % nv
    ! WRITE(200+THIS_PROC, *) N_ACC_ELEM+J,  &
    !                         ELEM(J) % V(1), ELEM(J) % V(2), ELEM(J) % V(3)
  end do

  ! Summ the buffer arrays for coordinates up
  call Global % Sum_Real_Array(nv_tot, Surf % buff_x)
  call Global % Sum_Real_Array(nv_tot, Surf % buff_y)
  call Global % Sum_Real_Array(nv_tot, Surf % buff_z)

  ! Summ the buffer arrays for elements' vertices up
  call Global % Sum_Int_Array(ne_tot, Surf % buff_i)
  call Global % Sum_Int_Array(ne_tot, Surf % buff_j)
  call Global % Sum_Int_Array(ne_tot, Surf % buff_k)
  call Global % Sum_Int_Array(ne_tot, Surf % buff_n)

  ! Fetch coordinates from all vertices
  do j = 1, nv_tot
    Vert(j) % x_n = Surf % buff_x(j)
    Vert(j) % y_n = Surf % buff_y(j)
    Vert(j) % z_n = Surf % buff_z(j)
    ! WRITE(300+THIS_PROC, *)  J, VERT(J) % X_N, VERT(J) % Y_N, VERT(J) % Z_N
  end do

  ! Fetch all elements' vertices
  do j = 1, ne_tot
    Elem(j) % v(1) = Surf % buff_i(j)
    Elem(j) % v(2) = Surf % buff_j(j)
    Elem(j) % v(3) = Surf % buff_k(j)
    Elem(j) % nv   = Surf % buff_n(j)
    ! WRITE(400+THIS_PROC, *) J, ELEM(J) % V(1), ELEM(J) % V(2), ELEM(J) % V(3)
  end do

  Surf % n_verts = nv_tot
  Surf % n_elems = ne_tot

  ! From this point on, all processors should have
  ! all the information on the entire grid

  end subroutine
