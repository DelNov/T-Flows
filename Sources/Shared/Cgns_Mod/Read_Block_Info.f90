!==============================================================================!
  subroutine Cgns_Mod_Read_Block_Info(base, block)
!------------------------------------------------------------------------------!
!   Gets n_bases from base node                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id             ! base index number    
  integer           :: block_id            ! block index number
  character(len=80) :: block_name          ! name of the block
  integer           :: block_mesh_info(3)  ! n_nodes, n_cells, and ...
                                           ! ... n_b_nodes(if sorted)
  integer           :: error
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block

  ! Read block information
  call Cg_Zone_Read_F(file_id,          &
                      base_id,          &
                      block_id,         &
                      block_name,       &
                      block_mesh_info,  &
                      error)

  if (error .ne. 0) then
    print *, '# Failed read block info'
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % name      = trim(block_name)
  cgns_base(base) % block(block) % mesh_info = block_mesh_info

  ! Total nodes and cells
  cnt_nodes = cnt_nodes + cgns_base(base) % block(block) % mesh_info(1)
  cnt_cells = cnt_cells + cgns_base(base) % block(block) % mesh_info(2)

  if(verbose) then
    print '(a)', ' #====================================================='
    print '(a,a36)', ' #     Block name: ',  &
             trim( cgns_base(base) % block(block) % name)
    print '(a)', ' #-----------------------------------------------------'
    print '(a,i35)', ' #     Block index: ', block
    print '(a,i41)', ' #     Nodes: ', &
      cgns_base(base) % block(block) % mesh_info(1)
    print '(a,i41)', ' #     Cells: ', &
      cgns_base(base) % block(block) % mesh_info(2)
    print '(a,i21)', ' #     Boundary nodes(if sorted): ',  &
      cgns_base(base) % block(block) % mesh_info(3)
  end if

  if (cgns_base(base) % block(block) % mesh_info(3) .ne. 0) then
    print *, '# Boundary condition nodes != 0 -> Unsupported'
    stop
  endif

  end subroutine
