!==============================================================================!
  subroutine Determine_Threads(Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is designed for configuring and managing threads in OpenMP
!>  for parallel processing within the Process. It determines the optimal
!>  distribution of computational workload across available threads based on
!>  the grid's structure and the specified thread count.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Setting Thread Count:                                                    !
!     - Determines the number of threads to be used based on user input and    !
!       system capabilities. Warns if the requested number of threads exceeds  !
!       the available resources.                                               !
!   * Cells and Faces Distribution:                                            !
!     - Distributes cells and faces across the available threads, aiming to    !
!       balance the computational load. This includes arranging cell threads   !
!       and aligning face threads with the corresponding cell threads.         !
!   * Data Reorganization for Parallel Processing:                             !
!     - Reorganizes internal data structures (such as cell and face indices)   !
!       to optimize them for parallel execution.                               !
!   * Debugging and Validation:                                                !
!     - Performs checks to ensure correct thread assignment and data           !
!       organization.                                                          !
!   * Performance Consideration:                                               !
!     - While testing with OpenMP has shown performance improvements, the      !
!       extent of these improvements may vary. The Incomplete Cholesky         !
!       preconditioner remains a significant factor in overall performance,    !
!       and users should manage their expectations regarding parallel          !
!       processing gains.  More work is needed in this direction.              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), target, intent(inout) :: Grid  !! The computational grid
!-----------------------------------[Locals]-----------------------------------!
  type(Vect_Type), pointer :: Vect
  integer                  :: c, s, m, i_fac
  integer                  :: n_cells_in, in_thread, n_remains, thr
  integer, allocatable     :: cells_in_thread(:)
  character(SL)            :: st1, st2
  integer                  :: f_s, l_s
  integer, allocatable     :: old_fc  (:,:)
  integer, allocatable     :: old_nn  (:)
  integer, allocatable     :: old_shad(:)
  integer, allocatable     :: old_nods(:,:)
!==============================================================================!

  m = size(Grid % faces_n, 1)
  allocate(old_fc  (2, Grid % n_faces))  ! old number of nodes
  allocate(old_nn  (   Grid % n_faces))  ! old number of nodes
  allocate(old_shad(   Grid % n_faces))
  allocate(old_nods(m, Grid % n_faces))

  ! Take aliases
  Vect => Grid % Vect

  !----------------------------!
  !                            !
  !   Find number of threads   !
  !                            !
  !- - - - - - - - - - - - - - +-----------------------!
  !   This is all a bit silly at this stage, it will   !
  !    make more sense when Vect_Mod matures a bit.    !
  !----------------------------------------------------!

  ! You read desired number of threads from control file, try to use that
  if(Vect % d_threads .lt. Vect % Get_Max_Threads()) then
    call Vect % Set_Num_Threads(Vect % d_threads)

  else if(Vect % d_threads .gt. Vect % Get_Max_Threads()) then
    write(st1, '(i0.0)') Vect % d_threads
    write(st2, '(i0.0)') Vect % Get_Max_Threads()
    call Message % Error(80, &
             'You are trying to run with '//trim(st1)//' OpenMP threads '  //  &
             'but there are only '//trim(st2)//' available. '              //  &
             'The workaround is to either increase the environment '       //  &
             'variable OMP_NUM_THREAD, or decrease parameter MAX_THREADS ' //  &
             'in control file.',                                               &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  ! Set n_threads to whatever you adjusted above
  Vect % n_threads = Vect % Get_Max_Threads()

  if(Vect % n_threads > 1) then
    write(st1, '(i0.0)') Vect % n_threads
    call Message % Framed(48, 'NOTE from OpenMP!',                  &
                              'You are running a simulation on ' // &
                               trim(st1) // ' threads.', one_proc=.true.)
  end if

  !-------------------------------------------------------!
  !                                                       !
  !   If there is more than one thread, distribute them   !
  !                                                       !
  !-------------------------------------------------------!
  if(PROGRAM_NAME == "Process" .and.  &
     Vect % n_threads > 1) then

    !--------------------------!
    !   Arrange cell threads   !
    !--------------------------!
    allocate(Vect % thread % f_cell(Vect % n_threads))
    allocate(Vect % thread % l_cell(Vect % n_threads))

    n_cells_in = Grid % region % l_cell(Grid % n_regions) -    &
                 Grid % region % f_cell(Grid % n_regions) + 1
    in_thread = n_cells_in / Vect % n_threads
    n_remains = mod(n_cells_in, Vect % n_threads)

    allocate(cells_in_thread(Vect % n_threads))
    cells_in_thread(:) = in_thread

    ! Add the remaining cells to the threads
    do thr = 1, n_remains
      cells_in_thread(thr) = cells_in_thread(thr) + 1
    end do

    ! Work out starts and ends of the threads
    Vect % thread % f_cell(1) = Grid % region % f_cell(Grid % n_regions)
    Vect % thread % l_cell(1) = Vect % thread % f_cell(1)  &
                              + cells_in_thread(1) - 1
    do thr = 2, Vect % n_threads
      Vect % thread % f_cell(thr) = Vect % thread % l_cell(thr-1) + 1
      Vect % thread % l_cell(thr) = Vect % thread % f_cell(thr)  &
                                  + cells_in_thread(thr) - 1
    end do

    ! Assign cell threads (not sure if really needed)
    do thr = 1, Vect % n_threads
      do c = Vect % thread % f_cell(thr), Vect % thread % l_cell(thr)
        Vect % cell_thread(c) = thr
      end do
    end do

    !-------------------------!
    !   Arange face threads   !
    !-------------------------!

    ! A plain copy of cell threads to face threads seems to work just fine
    ! (For the time being, I am only treating inside faces.  I am not even
    !  sure if explcit parallelization of boundary regions makes sense)
    do s = Faces_In_Domain_And_At_Buffers()
      c = Grid % faces_c(1, s)
      Vect % face_thread(s) = Vect % cell_thread(c)
    end do
    call Grid % Save_Vtu_Faces((/This_Proc(), N_Procs()/),  &
                               int_phi_f=Vect % face_thread)

    ! The problem with simple assignments above is that face threads
    ! are not continous, they break at buffer faces which are sorted
    ! by processors.  This shouldn't be needed for parallel runs ...
    do s = 1, Grid % n_faces
      Grid % old_f(s) = s
    end do

    ! ... hence, faces are sorted by threads
    f_s = Grid % region % f_face(Grid % n_regions)
    l_s = Grid % region % l_face(Grid % n_regions)
    call Sort % Int_Carry_Int(Vect % face_thread(f_s:l_s),  &
                                    Grid % old_f(f_s:l_s))

    ! Sort connectivity: faces_c, faces_n_nodes, faces_n and faces_s
    do s = Faces_In_Domain_And_At_Buffers()
      old_fc  (1:2,s) = Grid % faces_c  (1:2, Grid % old_f(s))
      old_nn  (    s) = Grid % faces_n_nodes( Grid % old_f(s))
      old_nods(1:m,s) = Grid % faces_n  (1:m, Grid % old_f(s))
      old_shad(    s) = Grid % faces_s      ( Grid % old_f(s))
    end do
    do s = Faces_In_Domain_And_At_Buffers()
      Grid % faces_c(1:2,  s) = old_fc  (1:2, s)
      Grid % faces_n_nodes(s) = old_nn  (     s)
      Grid % faces_n(1:m,  s) = old_nods(1:m, s)
      Grid % faces_s      (s) = old_shad(     s)
      Assert(Grid % faces_c(1,s) .lt. Grid % faces_c(2,s))
    end do

    ! Form the new face numbers from the old face numbers
    do s = 1, Grid % n_faces
      Grid % new_f(Grid % old_f(s)) = s
    end do

    ! Do the sorting of geometrical quantities. All ...
    ! ... which is read from .dim file should be here.
    call Sort % Real_By_Index(Grid % n_faces, Grid % sx(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % sy(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % sz(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % dx(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % dy(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % dz(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % f (1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % xf(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % yf(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % zf(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % rx(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % ry(1), Grid % new_f(1))
    call Sort % Real_By_Index(Grid % n_faces, Grid % rz(1), Grid % new_f(1))

    ! Correct shadow faces
    do s = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
      Grid % faces_s(s) = Grid % new_f(Grid % faces_s(s))
    end do

    ! Correct face indexes for cells
    do c = -Grid % n_bnd_cells, Grid % n_cells
      do i_fac = 1, Grid % cells_n_faces(c)
        s = Grid % cells_f(i_fac, c)             ! take the face's index
        Grid % cells_f(i_fac, c) = Grid % new_f(s)
      end do
    end do

    ! Work out starts and ends of the threads
    allocate(Vect % thread % f_face(Vect % n_threads))
    allocate(Vect % thread % l_face(Vect % n_threads))

    do thr = 1, Vect % n_threads
      Vect % thread % f_face(thr) =  HUGE_INT
      Vect % thread % l_face(thr) = -HUGE_INT
      do s = Faces_In_Domain_And_At_Buffers()
        if(Vect % face_thread(s) .eq. thr) then
          Vect % thread % f_face(thr) = min(Vect % thread % f_face(thr), s)
          Vect % thread % l_face(thr) = max(Vect % thread % l_face(thr), s)
        end if
      end do
    end do

    ! Perform some check
    do thr = 1, Vect % n_threads
      Assert(Vect % thread % f_face(thr) .le. Vect % thread % l_face(thr))
    end do
    do thr = 2, Vect % n_threads
      Assert(Vect % thread % f_face(thr) .ge. Vect % thread % f_face(thr-1))
      Assert(Vect % thread % l_face(thr) .ge. Vect % thread % l_face(thr-1))
      Assert(Vect % thread % f_face(thr) .ge. Vect % thread % l_face(thr-1))
    end do

  end if

  end subroutine
