!==============================================================================!
  subroutine Vof_Divide_Cuboid(Vof, c)
!------------------------------------------------------------------------------!
!   Divides cuboid recursively to make an accurate initialization              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: prelim_vof     => r_cell_01,  &
                      min_max_crit_1 => r_cell_02,  &
                      min_max_crit_2 => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type), target :: Vof
  integer                :: c
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: n, fu, n_node_vol, s
  integer                  :: ee, n_cylinders, iterat
  integer                  :: alloc_nod, alloc_fac
  integer                  :: i_nod, i_fac
  integer                  :: level_r
  real                     :: radius, height
  real                     :: p1_x, p1_y, p1_z
  real                     :: p2_x, p2_y, p2_z
  real                     :: res_dummy
  real                     :: seed_node(3)
  real,    allocatable     :: nod_x(:),nod_y(:),nod_z(:)
  integer, allocatable     :: faces_nod(:,:)
!==============================================================================!

  ! First take aliasesd
  Grid => Vof % pnt_grid
  n_nod_vol = 8
  alloc_nod = 0
  alloc_fac = 0

  level_r = 1
  do iterat = 0, level_r
    alloc_nod = alloc_nod + n_nod_vol * 8 ** iterat
  end do

  do iterat = 0, level_r
    alloc_fac = alloc_fac + 6 * 8 ** iterat
  end do

  allocate(nod_x(alloc_nod), nod_y(alloc_nod), nod_z(alloc_nod))
  nod_x = 0.0; nod_y = 0.0; nod_z = 0.0
  allocate(faces_nod(alloc_fac,4))

  !list nodes first cuboid:
  i_nod = 1
  do n = 1, Grid % cells_n_nodes(c)
    nod_x(inod) = Grid % xn(Grid % cells_n(n,c)) 
    nod_y(inod) = Grid % yn(Grid % cells_n(n,c))
    nod_z(inod) = Grid % zn(Grid % cells_n(n,c))
    i_nod = i_nod + 1
  end do

  !list of nodes on faces:
  i_fac = 1
  do iterat = 1, Grid % cells_n_faces(c)
    s = Grid % cells_f(iterat, c)
    do n = 1, Grid % faces_n_nodes(s)
      faces_nod(i_fac,n) = Grid % faces_n(n,s)
      i_fac = i_fac + 1
    end do
  end do


  ! find center (which is also nod on the slicing plane)
  seed_node = 0.0; 
  do n = i_nod, 1, -1 
    seed_node(1) = seed_node(1) + nod_x(n) / 8.0 
    seed_node(2) = seed_node(2) + nod_y(n) / 8.0 
    seed_node(3) = seed_node(3) + nod_z(n) / 8.0 
  end do

  end subroutine
