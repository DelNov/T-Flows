!==============================================================================!
  subroutine Cgns_Mod_Read_Base_Info(base)
!------------------------------------------------------------------------------!
!   Reads main info from base node base_id
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: base
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id    ! base index number
  character(len=80) :: base_name  ! name of the base
  integer           :: cell_dim   ! cell dimensions (3->volume cell)
  integer           :: phys_dim   ! number of coordinates to create vector
  integer           :: error
!==============================================================================!

  ! Set input parameters
  base_id = base

  ! Read CGNS base information
  call Cg_Base_Read_F(file_id,    &
                      base_id,    &
                      base_name,  &
                      cell_dim,   &
                      phys_dim,   &
                      error)

  if (error .ne. 0) then
    print *, '# Failed to get base info'
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % name     = base_name
  cgns_base(base) % cell_dim = cell_dim
  cgns_base(base) % phys_dim = phys_dim

  ! Print some info
  if(verbose) then
    print '(a)', ' #============================'
    print '(a,a9)', ' #   Base name:      ', trim(base_name)
    print '(a,i9)', ' #   Cell dimension: ', cell_dim
    print '(a,i9)', ' #   Phys dimension: ', phys_dim
    print '(a)', ' #----------------------------'
  end if

  end subroutine
