!==============================================================================!
  subroutine Cgns_Mod_Read_Interface_Info(base, block, interface)
!------------------------------------------------------------------------------!
!   Reads general interface info                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block, interface
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id              ! base index number
  integer           :: block_id             ! block index number
  integer           :: interface_id         ! interface index number
  character(len=80) :: interface_name       ! name of the interface
  character(len=80) :: interface_name_short ! shortened name of the interface
  integer           :: interface_type       ! interface type
  integer           :: interface_location   ! interface location
  integer           :: interface_ptset_type ! interface node/cell placement
  integer           :: interface_n_nodes    ! interface nodes or cells count
  character(len=80) :: donor_name           ! name of the donor
  integer           :: donor_block_type     ! donor block type
  integer           :: donor_data_type      ! donor data type
  integer           :: donor_ptset_type     ! donor node/cell placement
  integer           :: donor_n_nodes        ! donor nodes or cells count
  integer           :: error
!==============================================================================!

  ! Set input parameters
  base_id      = base
  block_id     = block
  interface_id = interface

  call Cg_Conn_Info_F(file_id,               & !(in )
                      base_id,               & !(in )
                      block_id,              & !(in )
                      interface_id,          & !(in )
                      interface_name,        & !(out)
                      interface_location,    & !(out)
                      interface_type,        & !(out)
                      interface_ptset_type,  & !(out)
                      interface_n_nodes,     & !(out)
                      donor_name,            & !(out)
                      donor_block_type,      & !(out)
                      donor_ptset_type,      & !(out)
                      donor_data_type,       & !(out)
                      donor_n_nodes,         & !(out)
                      error)                   !(out)

  if (error .ne. 0) then
    print *,"# Failed to read interface info"
    call Cg_Error_Exit_F()
    stop
  end if

  ! Cut substring from last "dom-"
  interface_name_short = interface_name(index(interface_name,"dom-", &
    back = .true.):)

  if(verbose) then
    print *, '#       ----------------------------------------'
    print '(2a)', ' #     Interface name: ', interface_name
    print '(2a)', ' #     Shortened to:   ', interface_name_short
    print '(2a)', ' #     Donor name:     ', donor_name
  end if

  ! Fetch received parameters
  cgns_base(base_id) % block(block_id) % interface(interface_id) % name = &
    trim(interface_name_short)

  end subroutine