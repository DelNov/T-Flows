!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Coordinates_In_Block(base, block)
!------------------------------------------------------------------------------!
!   Reads number of coordinates arrays from block_id
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer :: base_id         ! base index number
  integer :: block_id        ! block index number
  integer :: block_n_coords  ! number of coordinates in the block
  integer :: error           ! error status
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block

  ! Get number of coordinate arrays (1 for unstructure, 3 for structured)
  call Cg_Ncoords_F(file_id,         & !(in )
                    base_id,         & !(in )
                    block_id,        & !(in )
                    block_n_coords,  & !(out)
                    error)             !(out)

  if (error .ne. 0) then
    print *, "# Failed to get number of coordinate arrays"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % n_coords = block_n_coords

  end subroutine