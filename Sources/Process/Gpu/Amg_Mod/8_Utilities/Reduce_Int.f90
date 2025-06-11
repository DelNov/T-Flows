!==============================================================================!
  subroutine Reduce_Int(Amg, array, new_size)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)      :: Amg
  integer, allocatable :: array(:)
  integer              :: new_size
!-----------------------------------[locals]-----------------------------------!
  integer, allocatable :: temp(:)
  integer              :: old_size
!==============================================================================!

  ! If not allocated, just allocated and leave
  if(.not. allocated(array)) then
    allocate(array(new_size))
    array(1:new_size) = 0
    return
  end if

  ! If small enough, nothing to do
  old_size = size(array,1)
  if(old_size .le. new_size) return

  ! Store (a part of) the contents of the array
  allocate(temp(new_size))
  temp(1:new_size) = array(1:new_size)

  ! Resize the target array
  deallocate(array)
  allocate(array(new_size))

  ! Retreive the old data
  array(1:new_size) = temp(1:new_size)

  end subroutine
