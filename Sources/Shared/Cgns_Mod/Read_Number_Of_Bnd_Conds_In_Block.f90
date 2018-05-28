!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Bnd_Conds_In_Block(base, block)
!------------------------------------------------------------------------------!
!   Opens name_in file and return file index                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer :: base_id      ! base index number
  integer :: block_id     ! block index number
  integer :: n_bnd_conds  ! number of boundary conditions in block
  integer :: error
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block

  ! Get number of boundary condition in block
  call Cg_Nbocos_F(file_id,      &
                   base_id,      &
                   block_id,     &
                   n_bnd_conds,  &
                   error)   
  if (error.ne.0) then
    print *,"# Failed to obtain number of boundary conditions"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % n_bnd_conds = n_bnd_conds

  if(verbose) then
    print *, "#     Number of boundary conditions in the block:  ",  &
             cgns_base(base) % block(block) % n_bnd_conds
  end if
 
  allocate( cgns_base(base) % block(block) % bnd_cond(n_bnd_conds) )

  end subroutine
