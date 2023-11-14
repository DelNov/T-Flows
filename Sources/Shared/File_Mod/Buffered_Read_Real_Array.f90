!==============================================================================!
  subroutine Buffered_Read_Real_Array(File, file_unit, array)
!------------------------------------------------------------------------------!
!>  Reads real numbers from a file into an array using buffering. Buffering is
!>  beneficial when dealing with large data sets as it reduces the number of
!>  I/O operations by reading chunks of data at a time.
!------------------------------------------------------------------------------!
!   Functionality:
!
!   * Initialization:
!     - Sets tot to 0, and calculates lo, up, and n based on the bounds of
!       array.
!   * Reading Loop:
!     - The subroutine enters a loop to read the data in chunks (buffering).
!     - In each iteration, it determines rem, the number of elements to read
!       in this iteration, based on the remaining elements and the buffer size
!       (BUFFER_SIZE).
!     - It then reads rem number of real numbers from the file into
!       File % i_buffer.
!     - After reading into the buffer, it transfers the data from
!       File % i_buffer to array.
!   * Data Transfer:
!     - Uses a nested loop to transfer data from the buffer (File % i_buffer)
!       to array. The index j is calculated to correctly map buffer indices to
!       array indices.
!   * Check Completion:
!     - After each read, it updates tot and checks if all required data
!       (n elements) have been read. If so, it exits the loop.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)     :: File       !! parent class
  integer, intent(in)  :: file_unit  !! file unit
  real,    intent(out) :: array(:)   !! array to be read from the file
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n, tot, rem, lo, up
!==============================================================================!

  tot = 0
  lo  = lbound(array, 1)
  up  = ubound(array, 1)
  n   = up - lo + 1

  do
    ! Determine how many elements to read in this iteration
    rem = min(BUFFER_SIZE, n - tot)

    ! Read into the buffer
    read(file_unit) File % r_buffer(1:rem)

    ! Store the File % buffered data in your data structure
    do i = 1, rem
      j = lo - 1 + tot + i
      array(j) = File % r_buffer(i)
    end do

    tot = tot + rem

    ! Check if we've read all data
    if (tot .ge. n) exit

  end do

  end subroutine

