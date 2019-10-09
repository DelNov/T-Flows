!==============================================================================!
  subroutine Cgns_Mod_Write_Coordinate_Array(base, block, coord, grid)
!------------------------------------------------------------------------------!
!   Writes grid coordinates (RealDouble) [sequential vesion]                   !
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
  integer           :: coord_id        ! coord index number
  character(len=80) :: coord_name
  integer           :: i               ! lower range index
  integer           :: j               ! upper range index
  integer           :: error           ! error status
  real              :: coordinates(grid % n_nodes)
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  coord_id   = coord
  coord_name = cgns_base(base_id) % block(block_id) % coord_name(coord_id)

  i = 1
  j = cgns_base(base_id) % block(block_id) % mesh_info(1)

  ! Fetch received parameters
  select case (coord_id)
    case (1)
      coordinates = grid % xn(1:grid % n_nodes)
    case (2)
      coordinates = grid % yn(1:grid % n_nodes)
    case (3)
      coordinates = grid % zn(1:grid % n_nodes)
  end select

  ! Write grid coordinates
  call Cg_Coord_Write_F(file_id,      & !(in )
                        base_id,      & !(in )
                        block_id,     & !(in )
                        RealDouble,   & !(in )
                        coord_name,   & !(in )
                        coordinates,  & !(in )
                        coord_id,     & !(out)
                        error)          !(out)

  if (error .ne. 0) then
         print *, '#         Failed to write: ', trim(coord_name)
     call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose) then
    print *, '#         Coord array: ', coord_name
  end if
  if(verbose.and.coord_id.eq.1) then
    print *, '#         Number of nodes: ', j - i + 1
    print *, '#         First node:', i
    print *, '#         Last node: ', j
  end if

  end subroutine
