!==============================================================================!
  subroutine Cgns_Mod_Write_Solution_Info(base, block, solution)
!------------------------------------------------------------------------------!
!   Writes main info to solution node and sets its sol_id  [sequential vesion] !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block, solution
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id   ! base index number
  integer           :: block_id  ! block index number
  integer           :: sol_id    ! element section index
  character(len=80) :: sol_name  ! name of the Elements_t node
  integer           :: sol_type  ! types of elements in the section
  integer           :: error
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  sol_id   = solution

  sol_name = trim(cgns_base(base_id)%block(block_id)%solution(sol_id)%name)
  sol_type = cgns_base(base_id)%block(block_id)%solution(sol_id)%sol_type

  ! Create a FlowSolution_t node 
  call Cg_Sol_Write_F(file_id,   & !(in )
                      base_id,   & !(in )
                      block_id,  & !(in )
                      sol_name,  & !(in )
                      sol_type,  & !(in )
                      sol_id,    & !(out)
                      error)       !(out)

  if (error .ne. 0) then
    print *, '# Failed to write solution info'
    call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose ) then
    print *, '#         ---------------------------------'
    print *, '#         Solution name: ',   sol_name
    print *, '#         ---------------------------------'
    print *, '#         Solution idx:    ', sol_id
    print *, '#         Solution type:   ', GridLocationName(sol_type)
  end if


  end subroutine
