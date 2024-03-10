!==============================================================================!
  subroutine Test_002
!------------------------------------------------------------------------------!
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real, allocatable  :: a(:), b(:)
  integer, parameter :: N = 800*800*800
  integer, parameter :: N_STEPS = 1200  ! spend enough time on device
  integer            :: step
  real               :: dot, ts, te
!==============================================================================!

  print '(a)', ' #================================================='
  print '(a)', ' # TEST 2: Performing a vector vector dot product'
  print '(a)', ' #================================================='

  ! 800^3 => 7.629 GB on device => was OK
  print '(a,i12)', ' # The problem size is: ', N

  print '(a)', ' #----------------------------------------------------'
  print '(a)', ' # Be careful with memory usage.  If you exceed the'
  print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
  print '(a)', ' # card has the program will become memory bound no'
  print '(a)', ' # matter how you wrote it, and it may even crash.'
  print '(a)', ' #----------------------------------------------------'

  print '(a)', ' # Creating two vectors'
  allocate(a(N))
  allocate(b(N))

  a(:) = 1.0
  b(:) = 2.0

  ! Copy vectors to the device
  call Gpu % Vector_Real_Copy_To_Device(a)
  call Gpu % Vector_Real_Copy_To_Device(b)

  !-----------------------------------------------!
  !   Performing a fake time loop on the device   !
  !-----------------------------------------------!
  print '(a,i6,a)', ' # Performing a vector vector dot product ',  &
                    N_STEPS, ' times'

  call cpu_time(ts)
  do step = 1, N_STEPS
    call Linalg % Vec_D_Vec(N, dot, a, b)
  end do
  call cpu_time(te)

  ! Destroy data on the device, you don't need them anymore
  call Gpu % Vector_Real_Destroy_On_Device(a)
  call Gpu % Vector_Real_Destroy_On_Device(b)

  ! Print result
  print '(a,es12.3)', ' Dot product is: ', dot
  print '(a,es12.3)', ' Correct result: ', a(1) * b(1) * N

  print '(a,f12.3,a)', ' # Time elapsed for TEST 2: ', te-ts, ' [s]'

  end subroutine
