!==============================================================================!
  subroutine Enlarge_Real_By_Half(Amg, array, new_size)
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type)      :: Amg
  real,    allocatable :: array(:)
  integer              :: new_size
!-----------------------------------[Locals]-----------------------------------!
  real,    allocatable :: temp(:)
  integer              :: old_size
!==============================================================================!

  ! If not allocated, just allocated and leave
  if(.not. allocated(array)) then
    allocate(array(new_size + new_size / 2))
    array(1:new_size + new_size / 2) = 0.0
    return
  end if

  ! If large enough, nothing to do
  old_size = size(array,1)
  if(old_size .ge. new_size) return

  ! Store the contents of the array
  allocate(temp(old_size))
  temp(1:old_size) = array(1:old_size)

  ! Resize the target array
  deallocate(array)
  allocate(array(new_size + new_size / 2))

  ! Retreive the old data
  array(1:old_size) = temp(1:old_size)
  array(old_size+1:new_size + new_size / 2) = 0.0

  end subroutine
