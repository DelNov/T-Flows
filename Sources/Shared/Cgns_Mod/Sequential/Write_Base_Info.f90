!==============================================================================!
  subroutine Cgns_Mod_Write_Base_Info(base)
!------------------------------------------------------------------------------!
!   Writes main info to base node and sets its base_id  [sequential vesion]    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id    ! base index number
  character(len=80) :: base_name  ! name of the base
  integer           :: cell_dim   ! cell dimensions (3->volume cell)
  integer           :: phys_dim   ! number of coordinates to create vector
  integer           :: error
!==============================================================================!

  ! Set input parameters
  base_id   = base
  base_name = cgns_base(base_id) % name
  cell_dim  = cgns_base(base_id) % cell_dim
  phys_dim  = cgns_base(base_id) % phys_dim

  ! Read CGNS base information
  call Cg_Base_Write_F(file_id,   & !(in )
                       base_name, & !(in )
                       cell_dim,  & !(in )
                       phys_dim,  & !(in )
                       base_id,   & !(out)
                       error)       !(out)

  if (error .ne. 0) then
    print *, '# Failed to get base info'
    call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose) then
    print *, '#   ============================'
    print *, '#   Base name: ',      base_name
    print *, '#   Base id: ',        base_id
    print *, '#   Cell dimension: ', cell_dim
    print *, '#   Phys dimension: ', phys_dim
    print *, '#   ============================'
  end if

  end subroutine
