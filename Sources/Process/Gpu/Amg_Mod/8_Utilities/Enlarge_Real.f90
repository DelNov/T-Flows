!==============================================================================!
  subroutine Enlarge_Real(Amg, array, new_size)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)      :: Amg
  real,    allocatable :: array(:)
  integer              :: new_size
!-----------------------------------[locals]-----------------------------------!
  real,    allocatable :: temp(:)
  integer              :: old_size
!==============================================================================!

  ! If not allocated, just allocated and leave
  if(.not. allocated(array)) then
    allocate(array(new_size))
    array(1:new_size) = 0.0
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
  allocate(array(new_size))

  ! Retreive the old data
  array(1:old_size) = temp(1:old_size)
  array(old_size+1:new_size) = 0.0

  end subroutine
