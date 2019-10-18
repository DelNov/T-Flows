!==============================================================================!
  subroutine Multiphase_Mod_Vof_Initialization_Ellipsoid(mult)
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
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type),  pointer :: grid
  integer                   :: c, n
  integer                   :: ee, n_ellipses
  real                      :: radius_1, radius_2, radius_3
  real                      :: cent_x, cent_y, cent_z, res_dummy
!==============================================================================!

  ! First take aliases
  grid => mult % pnt_grid

  prelim_vof = 0.0

  ! Taking extrema:
  min_max_crit_1(:) =  HUGE
  min_max_crit_2(:) = -HUGE

  ! Open file to read Ellipsoid parameters:
  open(9, file = 'Ellipsoid-parameters')

  if (this_proc < 2) print *, '# Reading the file: ', 'Ellipsoid-parameters'

  read (9,*) n_ellipses

  do ee = 1, n_ellipses

    read (9,*) radius_1, radius_2, radius_3
    read (9,*) cent_x, cent_y, cent_z

    do c = 1, grid % n_cells
      ! for every node:
      do n = 1, grid % cells_n_nodes(c)
        res_dummy = sqrt(                                               &
              ((grid % xn(grid % cells_n(n,c))-cent_x)/radius_1)**2.0   &
            + ((grid % yn(grid % cells_n(n,c))-cent_y)/radius_2)**2.0   &
            + ((grid % zn(grid % cells_n(n,c))-cent_z)/radius_3)**2.0)

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

  close(9)

  end subroutine
