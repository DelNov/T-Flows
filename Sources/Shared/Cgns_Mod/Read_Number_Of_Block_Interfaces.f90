!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Block_Interfaces(base, block)
!------------------------------------------------------------------------------!
! Reads number of interfaces in block -> n_interfaces
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer :: base_id       ! base index number
  integer :: block_id      ! block index number
  integer :: number_of_interfaces
  integer :: error
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block

  ! Gets number of interfaces
  call Cg_Nconns_F(file_id,              & !(in )
                   base,                 & !(in )
                   block,                & !(in )
                   number_of_interfaces, & !(out)
                   error)                  !(out)

  if (error.ne.0) then
    print *, '# Failed to number of interfaces in block ', block
    call Cg_Error_Exit_F()
  endif

  if(verbose) then
    print '(a,i4)', &
      ' #     Interfaces in block: ', number_of_interfaces
  end if

  cgns_base(base) % block(block) % n_interfaces = number_of_interfaces

  ! Initialization of "interface" structure
  allocate( cgns_base(base) % block(block) % interface(number_of_interfaces) )

  cgns_base(base) % block(block) % interface(1:number_of_interfaces) % &
    marked_for_deletion = .false.

  cgns_base(base) % block(block) % interface(1:number_of_interfaces) % id = 0

  cgns_base(base) % block(block) % interface(1:number_of_interfaces) % type_c = ''

  end subroutine