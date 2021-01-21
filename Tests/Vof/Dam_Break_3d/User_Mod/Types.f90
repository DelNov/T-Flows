!==============================================================================!
!   Introduce new types to be used with User_Mod                               !
!==============================================================================!

  type probe
    integer, allocatable :: s_probe(:)
    real,    allocatable :: s_coor(:,:)
    real,    allocatable :: hcoor(:,:)
    real,    allocatable :: s_vof(:)
  end type probe

  ! Probes for height
  type(probe), target :: probes(4)

  ! Probes for pressure
  integer, allocatable :: nod_probe(:)
  real,    allocatable :: p_probe(:), p_dist(:)

