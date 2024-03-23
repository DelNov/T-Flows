!==============================================================================!
  subroutine Form_Momentum_Matrix(Process, Flow, Grid, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)         :: Process
  type(Field_Type),    target :: Flow
  type(Grid_Type), intent(in) :: Grid
  real,  optional, intent(in) :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Con_Type), pointer :: Mcon
  type(Sparse_Val_Type), pointer :: Mval
  real,      contiguous, pointer :: val(:), v_m(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:), row(:)
  real,      contiguous, pointer :: visc(:), dens(:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: m12, urf
# if T_FLOWS_DEBUG == 1
  real, allocatable :: work(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Momentum_Matrix')

  ! Take some aliases
  Mcon => Flow % Nat % C
  Mval => Flow % Nat % M
  val  => Flow % Nat % M % val  ! +
  dia  => Flow % Nat % C % dia  ! +
  pos  => Flow % Nat % C % pos  ! +
  fc   => Flow % Nat % C % fc   ! +
! row  => Flow % Nat % C % row
  v_m  => Flow % v_m
  dens => Flow % density
  visc => Flow % viscosity
  urf  =  Flow % u % urf
  nz   =  Flow % Nat % C % nonzeros

  Assert(urf > 0.0)

  !$acc kernels
  do i = 1, nz
    val(i) = 0.0
  end do
  !$acc end kernels

  !--------------------------------------------------!
  !   Compute neighbouring coefficients over cells   !
  !--------------------------------------------------!

  !$acc parallel loop independent
  do c1 = Cells_In_Domain()  ! that should be sufficient

    !$acc loop seq
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)

      ! Coefficients inside the domain
      if(c2 .gt. 0) then
        m12 = 0.5 * (visc(c1)+visc(c2)) * fc(s)
        if(c1 .lt. c2) then
          val(pos(1,s)) = -m12
          val(pos(2,s)) = -m12
        end if
        val(dia(c1)) = val(dia(c1)) + m12

      ! Coefficients at the boundaries
      else
        reg = Grid % region % at_cell(c2)
        if(Grid % region % type(reg) .eq. WALL .or.  &
           Grid % region % type(reg) .eq. INFLOW) then
          m12 = visc(c1) * fc(s)
          val(dia(c1)) = val(dia(c1)) + m12
        end if
      end if
    end do
    !$acc end loop

  end do
  !$acc end parallel

  !------------------------------------!
  !   Take care of the unsteady term   !
  !------------------------------------!
  if(present(dt)) then
    !$acc parallel loop independent
    do c = Cells_In_Domain()
      val(dia(c)) = val(dia(c)) + dens(c) * Grid % vol(c) / dt
    end do
    !$acc end parallel
  end if

  !--------------------------------------------------------------!
  !   Store volume divided by central coefficient for momentum   !
  !   and refresh its buffers before discretizing the pressure   !
  !--------------------------------------------------------------!
  !$acc parallel loop independent
  do c = Cells_In_Domain()
    v_m(c) = Grid % vol(c) / val(dia(c))
  end do
  !$acc end parallel

  call Grid % Exchange_Inside_Cells_Real(v_m)

  !-------------------------------------!
  !   Part 1 of the under-relaxation    !
  !   (Part 2 is in Compute_Momentum)   !
  !-------------------------------------!
  !$acc parallel loop independent
  do c = Cells_In_Domain()
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$acc end parallel

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = Cells_In_Domain()
    ! or: work(c) = val(dia(c))
    ! or: work(c) = Mcon % row(c+1) - Mcon % row(c)
    work(c) = v_m(c)
  end do
  call Grid % Exchange_Inside_Cells_Real(work)
  call Grid % Save_Debug_Vtu("v_m",              &
                             inside_name="v_m",  &
                             inside_cell=work)
# endif

  call Profiler % Stop('Form_Momentum_Matrix')

  end subroutine
