!==============================================================================!
  subroutine Cgns_Mod_Read_Coordinate_Info(base, block, coord)
!------------------------------------------------------------------------------!
!   Reads coord_name of coord_id
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: base, block, coord
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id          ! base index number
  integer           :: block_id         ! block index number
  integer           :: coord_id            
  integer           :: coord_data_type     
  character(len=80) :: coord_name
  integer           :: error            ! error status
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  coord_id = coord
 
  ! Get info about coordinate coord_id
  call Cg_Coord_Info_F(file_id,          &
                       base_id,          &
                       block_id,         &
                       coord_id,         &
                       coord_data_type,  &
                       coord_name,       &
                       error)             
  if (error .ne. 0) then
    print *, "# Failed to get info in for coord_id"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % coord_name(coord) = coord_name

  end subroutine
