!==============================================================================!
!   Compilation:                                                               !
!                                                                              !
!   gfortran -fdefault-real-8 -c -cpp ../Shared/Const_Mod.f90                  !
!   gfortran -fdefault-real-8 -c -cpp ../Shared/Comm_Mod_Seq.f90               !
!   gfortran -fdefault-real-8 -c -cpp ../Shared/Profiler_Mod.f90               !
!   gfortran -o prof -fdefault-real-8 -cpp Const_Mod.o Comm_Mod_Seq.o  \       !
!                                          Profiler_Mod.o Check_Profiler.f90   !
!------------------------------------------------------------------------------!

!==============================================================================!
  subroutine Long_Chunk()
!------------------------------------------------------------------------------!
!   This piece of code is just to do something ... well, something long.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: N = 1333
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, j, k
  real, dimension(N,N) :: a, b, c
!==============================================================================!

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

  end subroutine

!==============================================================================!
  subroutine Sub_3()
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!==============================================================================!

  call Profiler % Start('Sub_3')

  call Long_Chunk()               ! do your own long stuff

  call Profiler % Stop('Sub_3')

  end subroutine

!==============================================================================!
  subroutine Sub_2()
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  call Profiler % Start('Sub_2')

  call Long_Chunk()               ! do your own long stuff
  call Sub_3                      ! branch to another procedure

  call Profiler % Stop('Sub_2')

  end subroutine

!==============================================================================!
  subroutine Sub_1()
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  call Profiler % Start('Sub_1')

  call Long_Chunk()               ! do your own long stuff
  call Sub_2                      ! branch to another procedure)

  call Profiler % Stop('Sub_1')

  end subroutine

!==============================================================================!
  program Check_Profiler
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  call Profiler % Start('Main')

  call Long_Chunk()               ! do your own long stuff
  call Sub_1                      ! branch to another procedure

  call Profiler % Stop('Main')
  call Profiler % Statistics(1)

  end program
