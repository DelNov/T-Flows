!==============================================================================!
  subroutine Cgns_Mod_Read_Coordinate_Array(base, block, coord, grid)
!------------------------------------------------------------------------------!
!   Read grid coordinates (RealDouble)                                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, coord
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id         ! base index number
  integer           :: block_id        ! block index number
  character(len=80) :: coord_name      
  integer           :: i               ! lower range index
  integer           :: j               ! upper range index
  real, allocatable :: coordinates(:)  ! array of coordinate values
  integer           :: error           ! error status
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  coord_name = cgns_base(base) % block(block) % coord_name(coord)

  i = 1
  j = cgns_base(base) % block(block) % mesh_info(1)

  allocate(coordinates(i:j))

  ! Read grid x coordinates
  call Cg_Coord_Read_F(file_id,      &
                       base,         &
                       block,        &
                       coord_name,   &
                       RealDouble,   &
                       i,            &
                       j,            &
                       coordinates,  &
                       error)

  if (error.ne.0) then
    print *, "# Failed to read DoubleReal Coord", coord_name
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  select case (coord)
    case (1)
      i = cnt_nodes + 1
      j = cnt_nodes + j
      grid % xn(i:j) = coordinates(:)
    case (2)
      i = cnt_nodes + 1
      j = cnt_nodes + j
      grid % yn(i:j) = coordinates(:)
    case (3)
      i = cnt_nodes + 1
      j = cnt_nodes + j
      grid % zn(i:j) = coordinates(:)
  end select

  deallocate(coordinates)

  end subroutine
