!==============================================================================!
  subroutine Comm_Mod_Load_Maps(grid)
!------------------------------------------------------------------------------!
!   Fill the cell map for sequential version.                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s
!==============================================================================!

  !----------------------------------------------------------------------!
  !   For sequential run, no need to read the map, just set dimensions   !
  !----------------------------------------------------------------------!
  nc_s  = grid % n_cells
  nb_s  = grid % n_bnd_cells
  nc_t  = nc_s
  nb_t  = nb_s

  !-------------------------------------!
  !   Global cell numbers for T-Flows   !
  !-------------------------------------!
  allocate(grid % comm % cell_glo(-nb_s:nc_s)) 

  ! Fill up global cell numbers
  do c = -nb_t, nc_t
    grid % comm % cell_glo(c) = c
  end do
  
  !-----------------------------------------!
  !   Global cell numbers for MPI mapping   !
  !-----------------------------------------!
  allocate(grid % comm % cell_map    (nc_s))  
  allocate(grid % comm % bnd_cell_map(nb_s))

  ! -1 is to start from zero, as needed by MPI functions
  do c = 1, nc_t
    grid % comm % cell_map(c) = c - 1
  end do
  
  ! -1 is to start from zero, as needed by MPI functions
  do c = 1, nb_t
    grid % comm % bnd_cell_map(c) = c - 1
  end do

  end subroutine
