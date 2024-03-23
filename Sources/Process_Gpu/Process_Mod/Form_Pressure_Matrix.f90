!==============================================================================!
  subroutine Form_Pressure_Matrix(Process, Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Discreetized system of momentum conservation equations:                    !
!                                                                              !
!     [M]{u} = {b}   [kg m / s^2]   [N]                                        !
!                                                                              !
!   Dimensions of certain variables:                                           !
!                                                                              !
!     M               [kg / s]                                                 !
!     u, v, w         [m / s]                                                  !
!     b               [kg m / s^2]     [N]                                     !
!     p, pp           [kg / (m s^2)]   [N / m^2]                               !
!     p%x, p%y, p%z   [kg / (m^2 s^2)]                                         !
!     v_flux          [m^3 / s]                                                !
!------------------------------------------------------------------------------!
!   Discretized pressure-Poisson equation reads:                               !
!                                                                              !
!     [A] {pp} = {b}     [m^3 / s]                                             !
!                                                                              !
!   Dimensions of certain variables:                                           !
!                                                                              !
!     A               [m^4 s / kg]                                             !
!     pp              [kg / (m s^2)]                                           !
!     p%x, p%y, p%z   [kg / (m^2 s^2)]                                         !
!     b               [m^3 / s]                                                !
!------------------------------------------------------------------------------!
  class(Process_Type)         :: Process
  type(Field_Type),    target :: Flow
  type(Grid_Type), intent(in) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  type(Sparse_Con_Type), pointer :: Acon
  type(Sparse_Val_Type), pointer :: Aval, Mval
  real,      contiguous, pointer :: val(:), v_m(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  integer                        :: s, c1, c2, c, i_cel, i, nz
  real                           :: a12
  real, allocatable              :: work(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Pressure_Matrix')

  ! Take some aliases
  Acon => Flow % Nat % C
  Aval => Flow % Nat % A
  Mval => Flow % Nat % M
  val  => Flow % Nat % A % val  ! +
  v_m  => Flow % Nat % M % v_m  ! +
  dia  => Flow % Nat % C % dia  ! +
  pos  => Flow % Nat % C % pos  ! +
  fc   => Flow % Nat % C % fc
  nz   =  Flow % Nat % C % nonzeros

  !$acc kernels
  do i = 1, nz
    val(i) = 0.0
  end do
  !$acc end kernels

  !--------------------------------------------------------------------!
  !   This is cell based and will not create race conditions on GPUs   !
  !--------------------------------------------------------------------!

  !$acc parallel loop independent
  do c1 = Cells_In_Domain()  ! going through buffers should not be needed

    !$acc loop seq
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        a12 = fc(s) * 0.5 * (v_m(c1) + v_m(c2))
        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a12
        end if
        val(dia(c1)) = val(dia(c1)) + a12
      end if
    end do
    !$acc end loop

  end do
  !$acc end parallel

  ! De-singularize the system matrix ... just like this, ad-hoc
  !$acc parallel loop independent
  do c = Cells_In_Domain()
    val(dia(c)) = val(dia(c)) * (1.0 + MICRO)
  end do
  !$acc end parallel

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = 1, Grid % n_cells
    work(c) = val(dia(c))
  end do
  call Grid % Save_Debug_Vtu("a_diagonal",              &
                             inside_name="a_diagonal",  &
                             inside_cell=work)
# endif

  call Profiler % Stop('Form_Pressure_Matrix')

  end subroutine
