!==============================================================================!
  subroutine Vof_Initialization_Box(mult)
!------------------------------------------------------------------------------!
!   Initialize as vof = 1 all cells inside a rectangular box defined by eight  !
!   points, sorted as shown in the schematic below. It should work for any     !
!   orientation of the box and possibly, it can be generalized to any convex   !
!   polyedral                                                                  !
!------------------------------------------------------------------------------!
!                                                                              !
!     8-----------7                                                            !
!    /|          /|                                                            !
!   5-----------6 |                                                            !
!   | |         | |                                                            !
!   | |         | |                                                            !
!   | |         | |                                                            !
!   | 4 - - - - | 3                                                            !
!   |/          |/                                                             !
!   1-----------2                                                              !
!                                                                              !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: prelim_vof     => r_cell_01,  &
                      inside_c       => i_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type),  target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  real                      :: p_xyz(8,3)
  real,         allocatable :: p(:,:)
  integer                   :: c, n, fu
  integer                   :: ee, n_boxes, n_p, np_count, i_fac
  integer                   :: trios(6,3)
  real                      :: n_xyz(6,3), v1aux(3), v2aux(3)
  real                      :: res_dummy
  real                      :: dd(6)
  character(20)             :: p_name_aux
!==============================================================================!

  ! First take aliases
  grid => mult % pnt_grid

  prelim_vof = 0.0

  p_name_aux = 'obstacle'
  call Vof_Init_Random_Seed(p_name_aux)

  allocate(p(6,3))

  ! Initialize sorting for normals:
  trios(1,:) = (/1, 4, 2/)
  trios(2,:) = (/6, 7, 8/)
  trios(3,:) = (/3, 7, 6/)
  trios(4,:) = (/5, 8, 4/)
  trios(5,:) = (/1, 2, 6/)
  trios(6,:) = (/3, 4, 8/)

  ! Open file to read Plane parameters:
  call File_Mod_Open_File_For_Reading('box_parameters.ini', fu)

  read(fu,*) n_boxes

  np_count = 1

  do ee = 1, n_boxes

    do n_p = 1, 8
      read (fu,*) p_xyz(n_p,:)
    end do

    ! Define Planes of the box:

    do i_fac = 1, 6

      v1aux(:) = p_xyz(trios(i_fac,2),:) - p_xyz(trios(i_fac,1),:)
      v2aux(:) = p_xyz(trios(i_fac,3),:) - p_xyz(trios(i_fac,1),:)

      n_xyz(i_fac,1) =    v1aux(2) * v2aux(3) - v1aux(3) * v2aux(2)
      n_xyz(i_fac,2) = - (v1aux(1) * v2aux(3) - v1aux(3) * v2aux(1))
      n_xyz(i_fac,3) =    v1aux(1) * v2aux(2) - v1aux(2) * v2aux(1)

      dd(i_fac) = n_xyz(i_fac,1) * p_xyz(trios(i_fac,1),1)   &
                + n_xyz(i_fac,2) * p_xyz(trios(i_fac,1),2)   &
                + n_xyz(i_fac,3) * p_xyz(trios(i_fac,1),3)

      p(i_fac,:) = p_xyz(trios(i_fac,1), :)
    end do

    do c = 1, grid % n_cells
      inside_c(c) = Check_Inside_Box(mult, p, dd, n_xyz, c = c)
    end do

    ! Simply interpolate linearly
    do c = 1, grid % n_cells
      if (inside_c(c) == 1) then
        !prelim_vof(c) = 1.0 - (min_max_crit_2(c) - dd)  &
        !              / (min_max_crit_2(c)-min_max_crit_1(c))
        call Vof_Interface_Box(mult,               &
                               c,                  &
                               n_xyz,              &
                               dd,                 &
                               prelim_vof(c))

      else if (inside_c(c) == 2) then
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
