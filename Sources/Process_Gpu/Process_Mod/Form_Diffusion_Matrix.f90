!==============================================================================!
  subroutine Form_Diffusion_Matrix(Proc, Flow, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)         :: Proc
  type(Field_Type),    target :: Flow
  real,              optional :: dt                 !! time step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Sparse_Type), pointer :: M
  integer                    :: c, s, c1, c2, reg
  real                       :: visc, dens, m12
  real, allocatable          :: work(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Form_Diffusion_Matrix')

  ! Take some aliases
  Grid => Flow % pnt_grid
  M    => Flow % Nat % M
  dens =  Flow % density
  visc =  Flow % viscosity

  !---------------------------!
  !   Discretize the matrix   !
  !---------------------------!
  M % val(:) = 0.0

  !--------------------------------------------------!
  !   Compute neighbouring coefficients over faces   !
  !     This will create race conditions on GPUs     !
  !     but you know it and better do it on CPUs     !
  !--------------------------------------------------!

  ! Browse though all faces inside domain
  do s = Grid % n_bnd_cells + 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    Assert(c1 .gt. 0)
    Assert(c2 .gt. 0)

    ! Calculate coeficients for the momentum matrix
    ! Units: ............
    m12 = visc * Grid % s(s) / Grid % d(s)
    Assert(m12 .gt. 0.0)

    M % val(M % pos(1,s)) = -m12
    M % val(M % pos(2,s)) = -m12
    M % val(M % dia(c1))  = M % val(M % dia(c1)) + m12
    M % val(M % dia(c2))  = M % val(M % dia(c2)) + m12
  end do

  ! Increase central coefficient for walls
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c1 .gt. 0)
        Assert(c2 .lt. 0)

        m12 = visc * Grid % s(s) / Grid % d(s)
        Assert(m12 .gt. 0.0)

        M % val(M % dia(c1))  = M % val(M % dia(c1)) + m12
      end do
    end if
  end do

  !------------------------------------!
  !   Take care of the unsteady term   !
  !------------------------------------!
  if(present(dt)) then
    do c = 1, Grid % n_cells
      M % val(M % dia(c)) = M % val(M % dia(c)) + dens * Grid % vol(c) / dt
    end do
  end if

  !---------------------------------------------------------!
  !   Store volume divided by central coeff. for momentum   !
  !   That is a novelty here, it doesn't exist in T-Flows   !
  !---------------------------------------------------------!
  do c = 1, Grid % n_cells
    M % v_m(c) = Grid % vol(c) / M % val(M % dia(c))
  end do

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = 1, Grid % n_cells
    work(c) = M % val(M % dia(c))
    ! or: work(c) = M % row(c+1) - M % row(c) ?
  end do
  call Grid % Save_Debug_Vtu("m_diagonal",              &
                             inside_name="m_diagonal",  &
                             inside_cell=work)
# endif

  call Profiler % Stop('Form_Diffusion_Matrix')

  end subroutine
