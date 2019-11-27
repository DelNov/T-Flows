!==============================================================================!
  subroutine Vof_Initialization_Cylinder(mult)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: prelim_vof     => r_cell_01,  &
                      min_max_crit_1 => r_cell_02,  &
                      min_max_crit_2 => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  integer                   :: c, n, fu
  integer                   :: ee, n_cylinders
  real                      :: radius, height
  real                      :: p1_x, p1_y, p1_z
  real                      :: p2_x, p2_y, p2_z
  real                      :: res_dummy
!==============================================================================!

  ! First take aliasesd
  grid => mult % pnt_grid

  prelim_vof = 0.0

  ! Open file to read cylinder parameters:
  call File_Mod_Open_File_For_Reading('cylinder_parameters.dat', fu)

  read(fu,*) n_cylinders

  do ee = 1, n_cylinders

    ! Taking extrema
    min_max_crit_1(:) =  HUGE
    min_max_crit_2(:) = -HUGE

    read(fu,*) radius
    read(fu,*) p1_x, p1_y, p1_z
    read(fu,*) p2_x, p2_y, p2_z

    height = sqrt((p1_x-p2_x)**2+(p1_y-p2_y)**2+(p1_z-p2_z)**2)

    do c = 1, grid % n_cells
      ! for every node:
      do n = 1, grid % cells_n_nodes(c)

        res_dummy = (   &
           ( (p2_y-p1_y)*(p1_z-grid % zn(grid % cells_n(n,c)))        &
            -(p1_y-grid % yn(grid % cells_n(n,c)))*(p2_z-p1_z))**2    &
          +( (p2_x-p1_x)*(p1_z-grid % zn(grid % cells_n(n,c)))        &
            -(p1_x-grid % xn(grid % cells_n(n,c)))*(p2_z-p1_z))**2    &
          +( (p2_x-p1_x)*(p1_y-grid % yn(grid % cells_n(n,c)))        &
            -(p1_x-grid % xn(grid % cells_n(n,c)))*(p2_y-p1_y))**2 )  &
          / (radius*height) ** 2

        min_max_crit_1(c)= min(res_dummy, min_max_crit_1(c)) 
        min_max_crit_2(c)= max(res_dummy, min_max_crit_2(c)) 
      end do
    end do

    ! Simply interpolate linearly: 
    do c = 1, grid % n_cells
      if (min_max_crit_1(c) < 1.0 .and. min_max_crit_2(c) > 1.0) then
          prelim_vof(c) = (1.0 - min_max_crit_1(c))  &
                        / (min_max_crit_2(c)-min_max_crit_1(c))
      else if (min_max_crit_2(c) <= 1.0) then
          prelim_vof(c) = 1.0
      end if
    end do

    ! Precision  
    do c = 1, grid % n_cells
      mult % vof % n(c) = max(prelim_vof(c),mult % vof % n(c))
    end do

  end do

  close(fu)

  end subroutine
