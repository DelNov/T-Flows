!==============================================================================!
  subroutine Form_Diffusion_Matrix(Proc, Flow, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  real,           optional :: dt                 !! time step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),       pointer :: Grid
  type(Sparse_Con_Type), pointer :: Mcon
  type(Sparse_Val_Type), pointer :: Mval
  real, contiguous,      pointer :: visc(:), dens(:)
  integer                        :: c, s, c1, c2, reg
  real                           :: m12
# if T_FLOWS_DEBUG == 1
    real, allocatable :: work(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Form_Diffusion_Matrix')

  ! Take some aliases
  Grid => Flow % pnt_grid
  Mcon => Flow % Nat % C
  Mval => Flow % Nat % M
  dens => Flow % density
  visc => Flow % viscosity

  !---------------------------!
  !   Discretize the matrix   !
  !---------------------------!
  Mval % val(:) = 0.0

  !--------------------------------------------------!
  !   Compute neighbouring coefficients over faces   !
  !     This will create race conditions on GPUs     !
  !     but you know it and better do it on CPUs     !
  !--------------------------------------------------!

  ! Browse though all faces inside domain
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    Assert(c1 .gt. 0)
    Assert(c2 .gt. 0)

    ! Calculate coeficients for the momentum matrix
    ! Units: ............
    m12 = 0.5 * (visc(c1)+visc(c2)) * Grid % s(s) / Grid % d(s)
    Assert(m12 .gt. 0.0)

    Mval % val(Mcon % pos(1,s)) = -m12
    Mval % val(Mcon % pos(2,s)) = -m12
    Mval % val(Mcon % dia(c1))  = Mval % val(Mcon % dia(c1)) + m12
    Mval % val(Mcon % dia(c2))  = Mval % val(Mcon % dia(c2)) + m12
  end do

  ! Increase central coefficient for walls
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c1 .gt. 0)
        Assert(c2 .lt. 0)

        m12 = visc(c1) * Grid % s(s) / Grid % d(s)
        Assert(m12 .gt. 0.0)

        Mval % val(Mcon % dia(c1))  = Mval % val(Mcon % dia(c1)) + m12
      end do
    end if
  end do

  !------------------------------------!
  !   Take care of the unsteady term   !
  !------------------------------------!
  if(present(dt)) then
    do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
      Mval % val(Mcon % dia(c)) = Mval % val(Mcon % dia(c))  &
                                + dens(c) * Grid % vol(c) / dt
    end do
  end if

  !--------------------------------------------------------------!
  !   Store volume divided by central coefficient for momentum   !
  !   and refresh its buffers before discretizing the pressure   !
  !--------------------------------------------------------------!
  do c = 1, Grid % n_cells
    Mval % v_m(c) = Grid % vol(c) / Mval % val(Mcon % dia(c))
  end do
  call Grid % Exchange_Cells_Real(Mval % v_m)

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = 1, Grid % n_cells
    ! or: work(c) = Mval % val(Mcon % dia(c))
    work(c) = Mcon % row(c+1) - Mcon % row(c)
  end do
  call Grid % Exchange_Inside_Cells_Real(work)
  call Grid % Save_Debug_Vtu("m_diagonal",              &
                             inside_name="M_Diagonal",  &
                             inside_cell=work)
# endif

  call Profiler % Stop('Form_Diffusion_Matrix')

  end subroutine
