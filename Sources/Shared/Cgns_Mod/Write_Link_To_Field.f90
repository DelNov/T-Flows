!==============================================================================!
  subroutine Write_Link_To_Field(base, block, solution, input_name)
!------------------------------------------------------------------------------!
!   Writes links for fields in FlowSolution_t node in DB                       !
!   This allows to write WallDistance and CellDelta once then use links to it  !
!   [sequential and parallel version]                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: base, block, solution
  character(len=*) :: input_name
!-----------------------------------[Locals]-----------------------------------!
  integer       :: base_id        ! base index number
  integer       :: block_id       ! block index number
  integer       :: solution_id    ! solution index
  character(SL) :: base_name      ! name of the block
  character(SL) :: block_name     ! name of the block
  character(SL) :: solution_name  ! name of the Elements_t node
  character(SL) :: name_in_file   ! local var
  character(SL) :: node_name
  integer       :: error
!==============================================================================!

  ! Set input parameters
  base_id       = base
  base_name     = cgns_base(base_id) % name
  block_id      = block
  block_name    = cgns_base(base_id) % block(block_id) % name
  solution_id   = solution
  solution_name = trim(cgns_base(base_id) % block(block_id) % &
                  solution(solution_id) % name)

  !---------------------------------!
  !   Go to FlowSolution_t DB node  !
  !---------------------------------!

  call Cg_Goto_F(file_id,           & !(in )
                 base_id,           & !(in )
                 error,             & !(out)
                'Zone_t',           & !(in )
                 block_id,          & !(in )
                 'FlowSolution_t',  & !(in )
                 solution_id,       & !(in )
                 'end')

  if (error .ne. 0) then
    print *, '# Failed to navigate to: ', 'FlowSolution_t'
    call Cg_Error_Exit_F()
  endif

  !------------------!
  !   Write a link   !
  !------------------!

  node_name = trim(input_name)
  name_in_file = '/'//trim(base_name)//'/'&
    //trim(block_name)//'/'//trim(solution_name)//'/'//trim(node_name)

  ! Create a link in 'file_with_mesh' to 'name_in_file' and name it 'node_name'
  call Cg_Link_Write_F(trim(node_name),       & ! (in )
                       trim(file_with_mesh),  & ! (in )
                       trim(name_in_file),    & ! (in )
                       error)                   ! (out)

  if (error .ne. 0) then
    print *, '# Failed to create a link to ', trim(file_with_mesh)
    call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose .and. this_proc < 2) then
    print *, '# ', trim(node_name), ' linked to : ',  &
      trim(file_with_mesh), ':', trim(name_in_file)
  end if

  end subroutine
