!------------------------------------------------------------------------------!
!   Compilation:                                                               !
!                                                                              !
!   gfortran -fdefault-real-8 -c ../Shared/Const_Mod.f90                       !
!   gfortran -fdefault-real-8 -c ../Shared/Comm_Mod_Seq.f90                    !
!   gfortran -fdefault-real-8 -c ../Shared/Profiler_Mod.f90                    !
!   gfortran -o prof -fdefault-real-8  Const_Mod.o Comm_Mod.o Profiler_Mod.o \ !
!             Check_Profiler.f90                                               !
!------------------------------------------------------------------------------!

!==============================================================================!
  subroutine Sub_3()
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  integer, parameter :: N = 1333
!------------------------------------------------------------------------------!
  integer              :: i, j, k
  real, dimension(N,N) :: a, b, c
!==============================================================================!

  call Profiler % Start('Sub_3')

  do i = 1, N
    do j = 1, N
      a(i,j) = real(i*j)
      b(i,j) = sqrt(real(i*j))
      c(i,j) = 0.0
    end do
  end do

  do i = 1, N
    do j = 1, N
      do k = 1, N
        c(i,j) = c(i,j) + a (i,k) * b (k,j)
      end do
    end do
  end do

  call Profiler % Stop('Sub_3')

  end subroutine

!==============================================================================!
  subroutine Sub_2()
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  integer, parameter :: N = 1333
!------------------------------------------------------------------------------!
  integer              :: i, j, k
  real, dimension(N,N) :: a, b, c
!==============================================================================!

  call Profiler % Start('Sub_2')

  do i = 1, N
    do j = 1, N
      a(i,j) = real(i*j)
      b(i,j) = sqrt(real(i*j))
      c(i,j) = 0.0
    end do
  end do

  do i = 1, N
    do j = 1, N
      do k = 1, N
        c(i,j) = c(i,j) + a (i,k) * b (k,j)
      end do
    end do
  end do

  call Sub_3

  call Profiler % Stop('Sub_2')

  end subroutine

!==============================================================================!
  subroutine Sub_1()
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  integer, parameter :: N = 1333
!------------------------------------------------------------------------------!
  integer              :: i, j, k
  real, dimension(N,N) :: a, b, c
!==============================================================================!

  call Profiler % Start('Sub_1')

  do i = 1, N
    do j = 1, N
      a(i,j) = real(i*j)
      b(i,j) = sqrt(real(i*j))
      c(i,j) = 0.0
    end do
  end do

  do i = 1, N
    do j = 1, N
      do k = 1, N
        c(i,j) = c(i,j) + a (i,k) * b (k,j)
      end do
    end do
  end do

  call Sub_2

  call Profiler % Stop('Sub_1')

  end subroutine

!==============================================================================!
  program Long
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  integer, parameter :: N = 1333
!------------------------------------------------------------------------------!
  integer              :: i, j, k
  real, dimension(N,N) :: a, b, c
!==============================================================================!

  call Profiler % Start('Main')

  do i = 1, N
    do j = 1, N
      a(i,j) = real(i*j)
      b(i,j) = sqrt(real(i*j))
      c(i,j) = 0.0
    end do
  end do

  do i = 1, N
    do j = 1, N
      do k = 1, N
        c(i,j) = c(i,j) + a (i,k) * b (k,j)
      end do
    end do
  end do

  call Sub_1

  call Profiler % Stop('Main')
  call Profiler % Statistics(1)

  end program
