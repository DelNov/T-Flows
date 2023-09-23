!==============================================================================!
  subroutine Buffered_Read_Real_Array(File, file_unit, array)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)     :: File
  integer, intent(in)  :: file_unit
  real,    intent(out) :: array(:)
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

