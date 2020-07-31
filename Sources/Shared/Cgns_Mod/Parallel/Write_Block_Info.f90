!==============================================================================!
  subroutine Cgns_Mod_Write_Block_Info(base, block)
!------------------------------------------------------------------------------!
!   Writes zone/block info in base node [parallel vesion]                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer       :: base_id             ! base index number
  integer       :: block_id            ! block index number
  character(SL) :: block_name          ! name of the block
  integer       :: block_mesh_info(3)  ! n_nodes, n_cells, and ...
                                       ! ... n_b_nodes(if sorted)
  integer       :: error
!==============================================================================!

  ! Set input parameters
  base_id         = base
  block_id        = block
  block_name      = trim(cgns_base(base_id) % block(block_id) % name)
  block_mesh_info = cgns_base(base_id) % block(block_id) % mesh_info

  ! Create and/or write to a zone node
  call Cg_Zone_Write_F(file_id,         & !(in )
                       base_id,         & !(in )
                       block_name,      & !(in )
                       block_mesh_info, & !(in )
                       Unstructured,    & !(in )
                       block_id,        & !(out)
                       error)             !(out)

  if (error .ne. 0) then
    print *, '# Failed to write block info'
    call Cgp_Error_Exit_F()
  endif

  if(verbose .and. this_proc.lt.2) then
    print *, '#     =========================================='
    print *, '#     Block name:                ', block_name
    print *, '#     Block index:               ', block_id
    print *, '#     Nodes:                     ', block_mesh_info(1)
    print *, '#     Cells:                     ', block_mesh_info(2)
    print *, '#     Boundary nodes(if sorted): ', block_mesh_info(3)
    print *, '#     ------------------------------------------'
  end if

  if (block_mesh_info(3) .ne. 0) then
    print *, '# Boundary condition nodes != 0 -> Unsupported'
    stop
  endif

  end subroutine
