!==============================================================================!
  subroutine Vof_Initialization_Box(Vof)
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
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  real                         :: p_xyz(8,3)
  real,            allocatable :: p(:,:)
  integer                      :: c, n, i, fu
  integer                      :: ee, n_boxes, n_p, np_count, i_fac, seeds(1024)
  integer                      :: trios(6,3)
  real                         :: n_xyz(6,3), v1aux(3), v2aux(3)
  real                         :: res_dummy
  real                         :: dd(6)
  character(20)                :: p_name_aux
  real,    contiguous, pointer :: prelim_vof(:)
  integer, contiguous, pointer :: inside_c(:)
!==============================================================================!

  call Work % Connect_Real_Cell(prelim_vof)
  call Work % Connect_Int_Cell (inside_c)

  ! First take aliases
  Grid => Vof % pnt_grid

  prelim_vof = 0.0

  ! Initialize random seed
  p_name_aux = 'obstacle'
  call random_seed(size=n)
  seeds = (/(i, i=1,n)/)
  call random_seed(put=seeds(1:n))

  allocate(p(6,3))

  ! Initialize sorting for normals:
  trios(1,:) = (/1, 4, 2/)
  trios(2,:) = (/6, 7, 8/)
  trios(3,:) = (/3, 7, 6/)
  trios(4,:) = (/5, 8, 4/)
  trios(5,:) = (/1, 2, 6/)
  trios(6,:) = (/3, 4, 8/)

  ! Open file to read Plane parameters:
  call File % Open_For_Reading_Ascii('box_parameters.ini', fu)

  call File % Read_Line(fu)
  read(line % tokens(1), *) n_boxes

  np_count = 1

  do ee = 1, n_boxes

    do n_p = 1, 8
      call File % Read_Line(fu)
      read(line % tokens(1), *) p_xyz(n_p, 1)
      read(line % tokens(2), *) p_xyz(n_p, 2)
      read(line % tokens(3), *) p_xyz(n_p, 3)
    end do

    ! Define Planes of the box
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

    do c = Cells_In_Domain_And_Buffers()
      inside_c(c) = Check_Inside_Box(Vof, p, dd, n_xyz, c = c)
    end do

    ! Simply interpolate linearly
    do c = Cells_In_Domain_And_Buffers()
      if (inside_c(c) == 1) then
        !prelim_vof(c) = 1.0 - (min_max_crit_2(c) - dd)  &
        !              / (min_max_crit_2(c)-min_max_crit_1(c))
        call Vof_Interface_Box(Vof,               &
                               c,                  &
                               n_xyz,              &
                               dd,                 &
                               prelim_vof(c))
        prelim_vof(c) = 0.001

      else if (inside_c(c) == 2) then
        prelim_vof(c) = 1.0
      end if
    end do

    ! Precision
    do c = Cells_In_Domain_And_Buffers()
      Vof % fun % n(c) = max(prelim_vof(c),Vof % fun % n(c))
    end do

  end do

  close(fu)

  call Work % Disconnect_Real_Cell(prelim_vof)
  call Work % Disconnect_Int_Cell (inside_c)

  end subroutine
