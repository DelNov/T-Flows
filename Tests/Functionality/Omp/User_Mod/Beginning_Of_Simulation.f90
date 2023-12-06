!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: ITERS = 6000
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer    :: Grid
  real, contiguous, pointer    :: c_val1(:), c_val2(:)
  real, contiguous, pointer    :: f_val1(:), f_val2(:), f_val(:), sur(:)
  integer, contiguous, pointer :: faces_c(:,:)
  real                         :: tmp
  integer                      :: c, c1, c2, s, i, thr
  character(SL)                :: st
!==============================================================================!

  call Profiler % Stop('Main')

  call Message % Framed(60,                                      &
                        'Greetings ...', '... from '//__FILE__,  &
                        one_proc=.true.)

  !$omp parallel
  print *, 'Hello from ', omp_get_thread_num()+1
  !$omp end parallel

  call Work % Connect_Real_Cell(c_val1, c_val2)
  call Work % Connect_Real_Face(f_val1, f_val2, f_val)

  Grid    => Flow % pnt_Grid
  faces_c => Grid % faces_c
  sur     => Grid % s

  do c = 1, Grid % n_cells
    call random_number(tmp);  c_val1(c) = 10.0 * tmp
    call random_number(tmp);  c_val2(c) = 10.0 * tmp
  end do
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    f_val1(s) = c_val1(c1)
    f_val2(s) = c_val2(c2)
  end do

  !-----------------------------------!
  !   Test implicit loop over faces   !
  !-----------------------------------!
  call Profiler % Start('Implicit Loop Over Faces')
  do i = 1, ITERS

    !$omp parallel do         &
    !$omp private(s, c1, c2)  &
    !$omp shared (faces_c, sur, f_val, c_val1, c_val2)
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = faces_c(1, s)
      c2 = faces_c(2, s)
      f_val(s) = c_val1(c1) * sur(s)  &
               + c_val2(c2) * sur(s)
    end do
    !$omp end parallel do
  end do
  call Profiler % Stop('Implicit Loop Over Faces')

  !-----------------------------------!
  !   Test explicit loop over faces   !
  !   (Important: set thr private)    !
  !-----------------------------------!
  call Profiler % Start('Explicit Loop Over Faces')
  do i = 1, ITERS

    !$omp parallel                 &
    !$omp private(s, c1, c2, thr)  &
    !$omp shared (faces_c, sur, f_val, c_val1, c_val2)
    thr = omp_get_thread_num() + 1
    do s = Grid % Omp % thread % f_face(thr),  &
           Grid % Omp % thread % l_face(thr)
      c1 = faces_c(1, s)
      c2 = faces_c(2, s)
      f_val(s) = c_val1(c1) * sur(s)  &
               + c_val2(c2) * sur(s)
    end do
    !$omp end parallel
  end do
  call Profiler % Stop('Explicit Loop Over Faces')

  !------------------------------------!
  !   Test unwrapped loop over faces   !
  !------------------------------------!
  call Profiler % Start('Unwrapped Loop Over Faces')
  do i = 1, ITERS

    !$omp parallel do         &
    !$omp private(s, c1, c2)  &
    !$omp shared (faces_c, f_val1, f_val2, c_val1, c_val2)
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = faces_c(1, s)
      c2 = faces_c(2, s)
      f_val1(s) = c_val1(c1)
      f_val2(s) = c_val2(c2)
    end do
    !$omp end parallel do

    !$omp parallel do private(s) shared (sur, f_val, f_val1, f_val2)
    do s = Faces_In_Domain_And_At_Buffers()
      f_val(s) = f_val1(s) * sur(s)  &
               + f_val2(s) * sur(s)
    end do
    !$omp end parallel do

  end do
  call Profiler % Stop('Unwrapped Loop Over Faces')

  call Work % Disconnect_Real_Cell(c_val1, c_val2)
  call Work % Disconnect_Real_Face(f_val1, f_val2, f_val)

  call Profiler % Statistics(indent=34)

  call Global % End_Parallel
  stop

  end subroutine
