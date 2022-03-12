!==============================================================================!
  program Two_Point_Correlation
!------------------------------------------------------------------------------!
!   Program to extract two point correlations from T-Flows' monitoring files   !
!                                                                              !
!   The algorithm is taken from here:                                          !
!   https://nptel.ac.in/courses/112104118/lecture-32/32-4_correlation_fn.htm   !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  character(len= 80)        :: file_name = 'rb_conv-monit.001'
  character(len=300)        :: one_line
  integer                   :: k, n, dum
  integer, parameter        :: N_POINTS     =    40
  integer, parameter        :: N_TIME_STEPS = 12000
  real, allocatable, target :: u(:,:), v(:,:), w(:,:), p(:,:), t(:,:)
  real, allocatable         :: r(:), s1(:), s2(:)
  real, pointer             :: phi(:,:)
!==============================================================================!

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  allocate(u(N_TIME_STEPS, N_POINTS))
  allocate(v(N_TIME_STEPS, N_POINTS))
  allocate(w(N_TIME_STEPS, N_POINTS))
  allocate(p(N_TIME_STEPS, N_POINTS))
  allocate(t(N_TIME_STEPS, N_POINTS))
  allocate(r (N_POINTS))
  allocate(s1(N_POINTS))
  allocate(s2(N_POINTS))

  !-------------------------------------------!
  !   Chose the variable for which you want   !
  !   to extract the two point correlations   !
  !-------------------------------------------!
  phi => u

  !-------------------------------------------!
  !   Read all monitoring points from files   !
  !-------------------------------------------!
  do k = 1, N_POINTS
    write(file_name(15:17), '(i3.3)'), k

    open(9, file=file_name)

      ! Skip first line (which just gives coordinates)
      read(9,'(a300)', end=1) one_line

      do n = 1, N_TIME_STEPS
        read(9, *) dum, u(n,k), v(n,k), w(n,k), p(n,k), t(n,k)
      end do

    close(9)
  end do

  !-------------------------------------!
  !   Work out two point correlations   !
  !-------------------------------------!
  r (:) = 0.0
  s1(:) = 0.0
  s2(:) = 0.0

  ! Accumulate over time steps
  do n = 1, N_TIME_STEPS

    ! For all points 
    do k = 1, N_POINTS
      r (k) = r (k) + phi(n,1) * phi(n,k)  ! variable
      s1(k) = s1(k) + phi(n,1)**2          ! first sum
      s2(k) = s2(k) + phi(n,k)**2          ! second sum
    end do
  end do

  ! Time average
  r (:) = r (:) / N_TIME_STEPS
  s1(:) = s1(:) / N_TIME_STEPS
  s2(:) = s2(:) / N_TIME_STEPS

  ! Print correlation
  do k = 1, N_POINTS
    print *, k, r(k) / sqrt(s1(k)*s2(k))
  end do

  ! Error trap if end of file was reached
1 continue

  end program
