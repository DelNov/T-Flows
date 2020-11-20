include '../User_Mod/Vof_Initialization_Cylinder.f90'
include '../User_Mod/Vof_Initialization_Ellipsoid.f90'
include '../User_Mod/Vof_Initialization_Plane.f90'
include '../User_Mod/Vof_Initialization_Box.f90'
include '../User_Mod/Vof_Init_Random_Seed.f90'
include '../User_Mod/Vof_Interface_Cylinder.f90'
include '../User_Mod/Vof_Interface_Ellipsoid.f90'
include '../User_Mod/Vof_Interface_Plane.f90'
include '../User_Mod/Vof_Interface_Box.f90'
include '../User_Mod/Check_Inside_Cell.f90'
include '../User_Mod/Check_Inside_Box.f90'
include '../User_Mod/Convert_Problem_Name_To_Integer.f90'
include '../User_Mod/Vof_Area_Square_Circle.f90'
include '../User_Mod/Vof_Exact_Cylinder.f90'
include '../User_Mod/Vof_Quick_Sort.f90'
include '../User_Mod/Intersection_Line_Face.f90'
include '../User_Mod/Interpolate_From_Nodes.f90'
!==============================================================================!
  subroutine User_Mod_Initialize_Variables(flow, turb, mult, swarm, sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  type(Solver_Type),     target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  real,             pointer :: dt
  integer                   :: c, c1, c2, s, i_probe, c_inters, n, code, i_s
  integer                   :: i, j
  real                      :: fs, min_dist, dist, glo_dist, epsloc

  ! Probes for pressure
  integer, parameter        :: N_PROBE = 8
  real                      :: x_probe(N_PROBE) = (/0.8245,0.8245,0.8245,   &
                                                    0.8245,0.8035,0.7635,   &
                                                    0.7235,0.6835/)
  real                      :: y_probe(N_PROBE) = (/0.5000,0.5000,0.5000,   &
                                                    0.5000,0.5000,0.5000,   &
                                                    0.5000,0.5000/)
  real                      :: z_probe(N_PROBE) = (/0.0210,0.0610,0.1010,   &
                                                    0.1410,0.1610,0.1610,   &
                                                    0.1610,0.1610/)
  ! Probes for height
  integer, parameter        :: N_PROBE_F = 4
  real                      :: x_probe_f(N_PROBE_F) = (/0.4960,0.9920,1.4880,  &
                                                        2.6380/)
  real                      :: y_probe_f(N_PROBE_F) = (/0.5000,0.5000,0.5000,  &
                                                        0.5000/)
  real                      :: z_probe_f(N_PROBE_F) = (/0.0000,0.0000,0.0000,  &
                                                        0.0000/)
  real                      :: pab(2,3)
  real                      :: pint(3)
  integer,      allocatable :: s_probe(:), sp_point(:)
  real,         allocatable :: s_coor(:,:)
  integer,      allocatable :: s_temp(:)
  real,         allocatable :: sc_temp(:,:)
  logical                   :: inters
!==============================================================================!

  ! Take aliases
  grid  => flow % pnt_grid
  vof   => mult % vof
  dt    => flow % dt

  epsloc = epsilon(epsloc)

  ! Initialize the whole domain as 0.0
  do c = 1, grid % n_cells
    vof % n(c) = 0.0
  end do

  ! Under a Plane:
  ! call Vof_Initialization_Plane(mult)
  ! Ellipsoid:
  !call Vof_Initialization_Ellipsoid(mult)
  ! Cylinder:
  !call Vof_Initialization_Cylinder(mult)
  ! Box:
  call Vof_Initialization_Box(mult)

  call Grid_Mod_Exchange_Cells_Real(grid, vof % n)

  ! Old value
  vof % o(:) = vof % n(:)

  ! At faces

  ! At boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
      vof % n(c2) = vof % n(c1)
    else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
    else
      vof % n(c2) = vof % n(c1)
    end if
  end do

  ! Finding closest nodal points to probes
  ! looking only in the obstacle
  allocate(p_dist(N_PROBE), nod_probe(N_PROBE),    &
           p_probe(N_PROBE))
  p_dist = HUGE
  nod_probe = -1
  do i_probe = 1, N_PROBE
    min_dist = HUGE
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(Grid_Mod_Bnd_Cond_Name(grid,c2) .eq. 'STEP') then
        do n = 1, grid % cells_n_nodes(c1)
          dist = sqrt((grid % xn(grid % cells_n(n,c1))-x_probe(i_probe)) ** 2 &
                     +(grid % yn(grid % cells_n(n,c1))-y_probe(i_probe)) ** 2 &
                     +(grid % zn(grid % cells_n(n,c1))-z_probe(i_probe)) ** 2)
          if (dist<min_dist) then
            nod_probe(i_probe) = grid % cells_n(n,c1)
            p_dist(i_probe) = dist
            min_dist = dist
          end if
        end do
      end if
    end do
  end do

  do i_probe = 1, N_PROBE
    glo_dist = p_dist(i_probe)
    call Comm_Mod_Global_Min_Real(glo_dist)

    if(.not. Math_Mod_Approx_Real( glo_dist,    &
                                   p_dist(i_probe), epsloc)) then
      nod_probe(i_probe) = -1
      p_dist(i_probe) = HUGE
    end if

  end do

  ! Finding points at faces for fluid height
  ! Storing faces s that intersect the probe in each processor, and also
  ! the intersection point

  do i_probe = 1, N_PROBE_F
    pab(1,:) = (/x_probe_f(i_probe),          &
                 y_probe_f(i_probe),          &
                 z_probe_f(i_probe) - 1.0/)
    pab(2,:) = (/x_probe_f(i_probe),          &
                 y_probe_f(i_probe),          &
                 z_probe_f(i_probe) + 20.0/)  ! 20.0 > domain's height
    if (allocated(s_probe)) deallocate(s_probe, s_coor, sp_point)
    allocate(s_probe(1), s_coor(1,3), sp_point(N_PROBE_F))

    c_inters = 0
    do s = 1, grid % n_faces
      call Intersection_Line_Face(flow, mult, s, pab, pint, inters)

      if (inters) then
        if (c_inters == 0) then
          c_inters = c_inters + 1
          s_probe(c_inters) = s
          s_coor(c_inters, :) = pint
        else
          if (allocated(s_temp)) then
            deallocate(s_temp, sc_temp)
          end if
          allocate(s_temp(size(s_probe)), sc_temp(size(s_coor,1),3))
          s_temp = s_probe
          sc_temp = s_coor
          deallocate(s_probe, s_coor)
          allocate(s_probe(size(s_temp)+1), s_coor(size(sc_temp,1) + 1,3))
          s_probe(1:size(s_temp))    = s_temp(1:size(s_temp))
          s_coor (1:size(sc_temp,1),:) = sc_temp(1:size(sc_temp,1),:)
          c_inters = c_inters + 1
          s_probe(c_inters)  = s
          s_coor (c_inters,:) = pint
        end if
      end if
    end do
    if (c_inters > 0) then
      allocate(probes(i_probe) % s_probe(c_inters))
      allocate(probes(i_probe) % s_coor (c_inters, 3))
      allocate(probes(i_probe) % s_vof  (c_inters))
      probes(i_probe) % s_probe(1:c_inters)   = s_probe(1:c_inters)
      probes(i_probe) % s_coor (1:c_inters,:) = s_coor(1:c_inters,:)
      probes(i_probe) % s_vof  (1:c_inters)   = 0.0
    end if
  end do

  ! Sort probs by height
  do i_probe = 1, N_PROBE_F
    if (allocated(probes(i_probe) % s_probe)) then
      write(*,*) this_proc,i_probe, N_PROBE_F,c_inters
      call Vof_Quick_Sort(probes(i_probe) % s_probe(1:c_inters),   &
                          probes(i_probe) % s_coor(:,:), 3, 1,     &
                          size(probes(i_probe) % s_probe))
    end if
  end do

  end subroutine
