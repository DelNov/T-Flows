!==============================================================================!
  subroutine Cgns_Mod_Read_Block_Type(base, block)
!------------------------------------------------------------------------------!
!   Prints block type                                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer :: base_id             ! base index number    
  integer :: block_id            ! block index number
  integer :: block_type          ! type of the block 
  integer :: error
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block

  ! Read type of block
  call Cg_Zone_Type_F(file_id,     &
                      base_id,     &
                      block_id,    &
                      block_type,  &
                      error)
  if (error .ne. 0) then
    print *, "# Failed to get block type"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % type = block_type

  if(verbose) then
    print *, "#     Block type is ",  &
             ZoneTypeName(cgns_base(base) % block(block) % type)
  end if

  if (cgns_base(base) % block(block) % type .eq. STRUCTURED) then
    print *, "# Structured CGNS meshed are unsupported"
    stop
  endif

  end subroutine
