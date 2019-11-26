!==============================================================================!
  subroutine Multiphase_Mod_Vof_Initialization(mult, init_type)
!------------------------------------------------------------------------------!
!   Initialize Volume Fraction                                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: c_d => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  integer                       :: init_type
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  real,             pointer :: dt
  real,             pointer :: vof_f(:)
  integer                   :: c, c1, c2, s
  integer                   :: donor, accept
  real                      :: fs, dot_prod
  real                      :: alfa_u, alfa_d, alfa_a, alfa_d_til, alfa_cbc
  real                      :: alfa_uq, gamma_f, alfa_f_til, signo
  real                      :: beta_f, prodmag, ang, cod
!==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => flow % pnt_grid
  vof    => mult % vof
  dt     => flow % dt
  vof_f  => mult % vof_f

 ! Initialize the whole domain as 0.0
  do c = 1, grid % n_cells
    vof % n(c) = 0.0
  end do

  select case (init_type)
    case (1) ! Under a Plane:
      call Multiphase_Mod_Vof_Initialization_Plane(mult)
    case (2) ! Ellipsoid:
      call Multiphase_Mod_Vof_Initialization_Ellipsoid(mult)
    case (3) ! Cylinder:
      call Multiphase_Mod_Vof_Initialization_Cylinder(mult)
  end select

  call Comm_Mod_Exchange_Real(grid, vof % n)

  ! Old value
  vof % o(:) = vof % n(:)

  ! Initialize properties:
  do c = 1, grid % n_cells
    density(c)   = vof % n(c)         * phase_dens(1)      &
                 + (1.0 - vof % n(c)) * phase_dens(2)
    viscosity(c) = vof % n(c)         * phase_visc(1)      &
                 + (1.0 - vof % n(c)) * phase_visc(2)
  end do
  call Comm_Mod_Exchange_Real(grid, density)
  call Comm_Mod_Exchange_Real(grid, viscosity)

  ! Initialize volume fraction at faces:
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    vof_f(s) = 0.5 * (vof % n(c1) + vof % n(c2))
  end do

  ! Initialize density at faces:
  do s = 1, grid % n_faces
    dens_face(s) = vof_f(s)         * phase_dens(1)     &
                 + (1.0 - vof_f(s)) * phase_dens(2)
  end do

  end subroutine
