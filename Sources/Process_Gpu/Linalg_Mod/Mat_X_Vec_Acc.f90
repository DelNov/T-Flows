!==============================================================================!
  subroutine Mat_X_Vec_Acc(Lin, n, nz, c, a_val, a_col, a_row, b)
!------------------------------------------------------------------------------!
!>  This subroutine computes sparse-matrix vector multiplication on a device,
!>  without checking if variables are present on the device.
!------------------------------------------------------------------------------!
!   Notes:                                                                     !
!                                                                              !
!   * This subroutine used to have directives:                                 !
!     !$acc data present(c, a_val, a_col, a_row, b)                            !
!     ...                                                                      !
!     !$acc end data                                                           !
!     around the loop, but there was a problem with that.  The "end data"      !
!     would destroy data on the device which I don't want in an iterative      !
!     procedure, and "data present" couldn't hang here without "end data"      !
!                                                                              !
!   * The main loop here used to much more sophisticated, like the one you     !
!     would find in the textbooks or what my ChatGPT Assistant would write:    !
!                                                                              !
!     !$acc parallel loop                                                      !
!     do i = 1, n                                                              !
!       temp = 0.0                                                             !
!       !$acc loop reduction(+:temp)                                           !
!       do ij = a_row(i), a_row(i+1) - 1                                       !
!         j = a_col(ij)                                                        !
!         temp = temp + a_val(ij) * b(j)                                       !
!       end do                                                                 !
!       c(i) = temp                                                            !
!     end do                                                                   !
!                                                                              !
!     but this version led to no improvement in performance whatsoever.        !
!     Once I started to use the simple:                                        !
!                                                                              !
!     !$acc kernels                                                            !
!     ...                                                                      !
!     !$acc end kernels                                                        !
!                                                                              !
!     speed up was there like in the case of full matrices.                    !
!                                                                              !
!   * Using intent clause here, was causing slower runs.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin         !! parent class
  integer            :: n           !! matrix and vector dimension
  integer            :: nz          !! number of nonzeros
  real               :: c(n)        !! result vector
  real               :: a_val(nz)   !! operand matrix values
  integer            :: a_col(nz)   !! operand matrix columns
  integer            :: a_row(n+1)  !! operand matrix rows
  real               :: b(n)        !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, ij
  real    :: temp
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Lin)
!==============================================================================!

  !$acc kernels
  do i = 1, n
    temp = 0.0
    do ij = a_row(i), a_row(i+1) - 1
      j = a_col(ij)
      temp = temp + a_val(ij) * b(j)
    end do
    c(i) = temp
  end do
  !$acc end kernels

  end subroutine

