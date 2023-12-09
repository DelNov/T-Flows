!==============================================================================!
  subroutine Allocate_Nodes(Grid, nn)
!------------------------------------------------------------------------------!
!>  Allocates memory for node-based data (arrays and matrices), for geometrical
!>  (xn, yn, zn ...) and connectivity data (new_n, old_n, Comm % node_glo).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! grid under consideration
  integer, intent(in) :: nn    !! number of nodes in the grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: n
!==============================================================================!

  ! Store the number of nodes for the grid
  Grid % n_nodes = nn

  ! Allocate memory for node coordinates
  allocate(Grid % xn(1:nn));  Grid % xn(:) = 0.0
  allocate(Grid % yn(1:nn));  Grid % yn(:) = 0.0
  allocate(Grid % zn(1:nn));  Grid % zn(:) = 0.0

  allocate(Grid % Comm % node_glo(1:nn));  Grid % Comm % node_glo(:) = 0
  do n = 1, nn
    Grid % Comm % node_glo(n) = n
  end do

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  if(PROGRAM_NAME .ne. "Process") then
    allocate(Grid % new_n(n));  Grid % new_n(:) = 0
    allocate(Grid % old_n(n));  Grid % old_n(:) = 0
  end if

  ! This is used only in initial stages of Generate, inside Domain_Mod really
  if(PROGRAM_NAME .eq. 'Generate') then
    allocate(Grid % node_positioned(n));  Grid % node_positioned(:) = .false.
  end if

  end subroutine
