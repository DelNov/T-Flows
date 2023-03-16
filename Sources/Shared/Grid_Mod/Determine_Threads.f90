!==============================================================================!
  subroutine Determine_Threads(Grid)
!------------------------------------------------------------------------------!
!   Take care of threads for OpenMP                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), target, intent(inout) :: Grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Vect_Type), pointer :: Vect
  integer                  :: c, s
  integer                  :: n_cells_in, in_region, n_remains, reg
  integer, allocatable     :: cells_in_region(:)
  character(SL)            :: st1, st2
!==============================================================================!

  ! Take aliases
  Vect => Grid % Vect

  !----------------------------!
  !   Find number of threads   !
  !- - - - - - - - - - - - - - +-----------------------!
  !   This is all a bit silly at this stage, it will   !
  !    make more sense when Vect_Mod matures a bit.    !
  !----------------------------------------------------!

  ! You read desired number of threads from control file, try to use that
  if(Grid % Vect % d_threads .lt. Grid % Vect % Get_Max_Threads()) then
    call Grid % Vect % Set_Num_Threads(Grid % Vect % d_threads)

  else if(Grid % Vect % d_threads .gt. Grid % Vect % Get_Max_Threads()) then
    write(st1, '(i0.0)') Grid % Vect % d_threads
    write(st2, '(i0.0)') Grid % Vect % Get_Max_Threads()
    call Message % Error(80, &
             'You are trying to run with '//trim(st1)//' OpenMP threads '  //  &
             'but there are only '//trim(st2)//' available. '              //  &
             'The workaround is to either increase the environment '       //  &
             'variable OMP_NUM_THREAD, or decrease parameter MAX_THREADS ' //  &
             'in control file.',                                               &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  ! Set n_threads to whatever you adjusted above
  Grid % Vect % n_threads = Grid % Vect % Get_Max_Threads()

  if(Grid % Vect % n_threads > 1) then
    write(st1, '(i0.0)') Grid % Vect % n_threads
    call Message % Framed(48, 'NOTE from OpenMP!',                  &
                              'You are running a simulation on ' // &
                               trim(st1) // ' threads.', one_proc=.true.)
  end if

  !----------------------------------------------------------!
  !   If there is more than one thread, decompose the mesh   !
  !----------------------------------------------------------!
  if(PROGRAM_NAME == "Process" .and.  &
     Grid % Vect % n_threads > 1) then

    allocate(Grid % Vect % region % f_cell(Grid % Vect % n_threads))
    allocate(Grid % Vect % region % l_cell(Grid % Vect % n_threads))

    n_cells_in = Grid % region % l_cell(Grid % n_regions) -    &
                 Grid % region % f_cell(Grid % n_regions) + 1

    in_region = n_cells_in / Grid % Vect % n_threads
    n_remains = mod(n_cells_in, Grid % Vect % n_threads)

    allocate(cells_in_region(Grid % Vect % n_threads))
    cells_in_region(:) = in_region

    ! Add the remaining cells to the regions
    do reg = 1, n_remains
      cells_in_region(reg) = cells_in_region(reg) + 1
    end do

    ! Work out starts and ends of the regions
    Vect % region % f_cell(1) = Grid % region % f_cell(Grid % n_regions)
    Vect % region % l_cell(1) = Vect % region % f_cell(1)  &
                              + cells_in_region(1) - 1
    do reg = 2, Grid % Vect % n_threads
      Vect % region % f_cell(reg) = Vect % region % l_cell(reg-1) + 1
      Vect % region % l_cell(reg) = Vect % region % f_cell(reg)  &
                                + cells_in_region(reg) - 1
    end do

    ! Assign cell threads (not sure if really needed)
    do reg = 1, Vect % n_threads
      do c = Vect % region % f_cell(reg), Vect % region % l_cell(reg)
        Vect % cell_thread(c) = reg
      end do
    end do

    ! Now this was good, but was having a hell of an impact
    ! on communication patterns and backup files creation
    ! call Grid % Sort_Cells_By_Thread()
    call Grid % Save_Vtu_Faces()

  end if

  end subroutine
