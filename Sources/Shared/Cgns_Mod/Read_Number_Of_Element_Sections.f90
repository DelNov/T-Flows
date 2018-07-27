!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Element_Sections(base, block)
!------------------------------------------------------------------------------!
!   Gets n_sects from block
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer :: base_id   ! base index number
  integer :: block_id  ! block index number
  integer :: n_sects   ! number of element sections in a block
  integer :: error
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block

  ! Get number of element sections
  call Cg_Nsections_F(file_id,   & !(in )
                      base_id,   & !(in )
                      block_id,  & !(in )
                      n_sects,   & !(out)
                      error)       !(out)

  if (error.ne.0) then
    print *, "# Failed to read number of elements"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % n_sects = n_sects

  if(verbose) then
    print *, "#       Number of sections: ",  &
             cgns_base(base) % block(block) % n_sects
  end if

  allocate( cgns_base(base) % block(block) % section(n_sects) )

  end subroutine