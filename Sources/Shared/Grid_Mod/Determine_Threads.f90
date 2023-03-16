!==============================================================================!
  subroutine Determine_Threads(Grid)
!------------------------------------------------------------------------------!
!   Take care of threads for OpenMP                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, s
  integer       :: c_f, c_l, s_f, s_l      ! first/last cell, first/last face
  character(SL) :: st1, st2
!==============================================================================!


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

    ! Set the face range to the inside faces
    s_f = Grid % region % f_face(Grid % n_regions)
    s_l = Grid % region % l_face(Grid % n_regions)

    ! In case of parallel runs, you do not want the faces in buffers,
    ! you want threads only in inside cells, not spreading to buffers
    do s = Faces_In_Domain()
      if(Grid % Comm % cell_proc(Grid % faces_c(2, s)) .ne. this_proc) then
        s_l = s - 1
        exit
      end if
    end do

    ! Check 1
    do s = s_f, s_l
      Assert(Grid % Comm % cell_proc(Grid % faces_c(1, s)) .eq. this_proc)
      Assert(Grid % Comm % cell_proc(Grid % faces_c(2, s)) .eq. this_proc)
    end do

    ! Create METIS
    call Metis % Create_Metis(s_f, s_l, Grid % faces_c, Grid % Vect % n_threads)

    ! Set some high value for threads everywhere.  In calls which follow
    ! cells which are not in buffers will be over-written by true threads
    Grid % Vect % cell_thread(:) = Grid % Vect % n_threads * 10

    ! Set the cell range to the inside of the domain, no buffer cells here
    c_f = Grid % region % f_cell(Grid % n_regions)
    c_l = Grid % region % l_cell(Grid % n_regions)

    ! Call METIS to assign regions to inside cells (not in buffers).  In
    ! buffers, you should still have values you assigned a few lines above
    call Metis % Call_Metis(Grid % Vect % n_threads,  &
                            Grid % Vect % cell_thread(c_f:c_l))

    do c = Cells_In_Domain()
      Assert(Grid % Vect % cell_thread(c) .le. Grid % Vect % n_threads)
    end do
    do c = Cells_In_Buffers()
      Assert(Grid % Vect % cell_thread(c) .gt. Grid % Vect % n_threads)
    end do

    ! Now this is quite brave
    call Grid % Sort_Cells_By_Thread()

  end if

  end subroutine
