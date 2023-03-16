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
  integer                  :: c, s, c1, c2
  integer                  :: n_cells_in, in_thread, n_remains, thr
  integer, allocatable     :: cells_in_thread(:)
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

  !----------------------------------------------------------!
  !   If there is more than one thread, decompose the mesh   !
  !----------------------------------------------------------!
  if(PROGRAM_NAME == "Process" .and.  &
     Vect % n_threads > 1) then

    allocate(Vect % thread % f_cell(Vect % n_threads))
    allocate(Vect % thread % l_cell(Vect % n_threads))

    n_cells_in = Grid % region % l_cell(Grid % n_regions) -    &
                 Grid % region % f_cell(Grid % n_regions) + 1

    in_thread = n_cells_in / Vect % n_threads
    n_remains = mod(n_cells_in, Vect % n_threads)

    allocate(cells_in_thread(Vect % n_threads))
    cells_in_thread(:) = in_thread

    ! Add the remaining cells to the regions
    do thr = 1, n_remains
      cells_in_thread(thr) = cells_in_thread(thr) + 1
    end do

    ! Work out starts and ends of the regions
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

    ! A plain copy of cell threads to face threads seems to work just fine
    ! (For the time being, I am only treating inside faces.  I am not even
    !  sure if explcit parallelization of boundary regions makes sense)
    do s = Faces_In_Domain()
      c = Grid % faces_c(1, s)
      Vect % face_thread(s) = Vect % cell_thread(c)
    end do
    call Grid % Save_Vtu_Faces(int_phi_f=Vect % face_thread)

    ! Now this was good, but was having a hell of an impact
    ! on communication patterns and backup files creation
    ! call Grid % Sort_Cells_By_Thread()

  end if

  end subroutine
