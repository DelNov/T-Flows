!==============================================================================!
  subroutine Test_003
!------------------------------------------------------------------------------!
!>  Tests vector operations of the form c = a +/- s * b
!------------------------------------------------------------------------------!
  use Gpu_Mod
  use Linalg_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real, allocatable  :: a(:), b(:), c(:), d(:)
  integer, parameter :: N = 600*600*600
  integer, parameter :: N_STEPS = 1200  ! spend enough time on device
  integer            :: step
!==============================================================================!

  ! Start the parallel run and the profiler
  call Global % Start_Parallel
  call Profiler % Start('Test_003')

  ! Check if it was run in parallel
  if(Parallel_Run()) then
    call Message % Error(60, 'This test make no sense in parallel, '//  &
                             'run it sequentially', one_proc=.true.)
    call Global % Wait()
    call Global % End_Parallel
    stop
  end if

  print '(a)', ' #===================================================='
  print '(a)', ' # TEST 3: Performing vector operations:'
  print '(a)', ' #         c = a + s * b  and  c = a - s * b'
  print '(a)', ' #===================================================='

  print '(a,i12)', ' # The problem size is: ', N

# if T_FLOWS_GPU == 1
    print '(a)', ' #----------------------------------------------------'
    print '(a)', ' # Be careful with memory usage.  If you exceed the'
    print '(a)', ' # 90% (as a rule of thumb) of the memory your GPU'
    print '(a)', ' # card has the program will become memory bound no'
    print '(a)', ' # matter how you wrote it, and it may even crash.'
    print '(a)', ' #----------------------------------------------------'
# endif

  print '(a)', ' # Creating four vectors'
  allocate(a(N))
  allocate(b(N))
  allocate(c(N))
  allocate(d(N))

  a(:) = 1.0
  b(:) = 2.0

  ! Copy vectors and create on the device
  call Gpu % Vector_Real_Copy_To_Device(a)
  call Gpu % Vector_Real_Copy_To_Device(b)
  call Gpu % Vector_Real_Create_On_Device(c)
  call Gpu % Vector_Real_Create_On_Device(d)

  !-----------------------------------------------!
  !   Performing a fake time loop on the device   !
  !-----------------------------------------------!
  print '(a,i6,a)', ' # Performing a sparse-matrix vector product',  &
                    N_STEPS, ' times'

  call Profiler % Start('Useful_Work')
  do step = 1, N_STEPS
    call Linalg % Vec_P_Sca_X_Vec(N, c, a,  2.0, b)  ! result should be  5
    call Linalg % Vec_P_Sca_X_Vec(N, d, c, -2.0, b)  ! result should be  1
  end do
  call Profiler % Stop('Useful_Work')

  ! Copy results back to host
  call Gpu % Vector_Update_Host(c)
  call Gpu % Vector_Update_Host(d)

  ! Destroy data on the device, you don't need them anymore
  call Gpu % Vector_Real_Destroy_On_Device(a)
  call Gpu % Vector_Real_Destroy_On_Device(b)
  call Gpu % Vector_Real_Destroy_On_Device(c)
  call Gpu % Vector_Real_Destroy_On_Device(d)

  ! Print results
  print '(a,es12.3)', ' vector c(1  ):', c(1  )
  print '(a,es12.3)', ' vector c(2  ):', c(2  )
  print '(a,es12.3)', ' vector c(N-1):', c(N-1)
  print '(a,es12.3)', ' vector c(N  ):', c(N  )
  print '(a,es12.3)', ' vector c(1  ):', d(1  )
  print '(a,es12.3)', ' vector c(2  ):', d(2  )
  print '(a,es12.3)', ' vector c(N-1):', d(N-1)
  print '(a,es12.3)', ' vector c(N  ):', d(N  )

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_003')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine
