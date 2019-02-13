!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Blocks_In_Base(base)
!------------------------------------------------------------------------------!
!   Gets n_blocks from base node base_id
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base
!-----------------------------------[Locals]-----------------------------------!
  integer :: base_id   ! base index number     
  integer :: n_blocks  ! number of blocks present in base
  integer :: error     ! error status
!==============================================================================!

  ! Set input parameters
  base_id = base

 ! Get number of blocks in base
  call Cg_Nzones_F(file_id,   & !(in )
                   base_id,   & !(in )
                   n_blocks,  & !(out)
                   error)       !(out)

  if (error .ne. 0) then
    print *, "# Failed to get blocks number"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % n_blocks = n_blocks

  ! Print some info
  if(verbose) then
    print "(a,i8)", " #   Number of blocks:", cgns_base(base) % n_blocks
  end if

  ! Allocate memory for the blocks in current base
  allocate(cgns_base(base) % block(cgns_base(base) % n_blocks))

  end subroutine
