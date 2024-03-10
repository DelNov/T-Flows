!==============================================================================!
  subroutine Form_Pressure_Matrix(Proc, Flow)
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
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Sparse_Type), pointer :: A, M
  integer                    :: s, c1, c2, c
  real                       :: a12
  real, allocatable          :: work(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Form_Pressure_Matrix')

  ! Take some aliases
  Grid => Flow % pnt_grid
  A    => Flow % Nat % A
  M    => Flow % Nat % M

  !----------------------------------------------!
  !   This will create race conditions on GPUs   !
  !----------------------------------------------!
  do s = Grid % n_bnd_cells + 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Calculate coeficients for the pressure matrix
    ! Units: m * m^3 s / kg = m^4 s / kg
    a12 = A % fc(s) * 0.5 * (M % v_m(c1) + M % v_m(c2))
    A % val(A % pos(1,s)) = -a12
    A % val(A % pos(2,s)) = -a12
    A % val(A % dia(c1))  = A % val(A % dia(c1)) + a12
    A % val(A % dia(c2))  = A % val(A % dia(c2)) + a12

  end do

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = 1, Grid % n_cells
    work(c) = A % val(A % dia(c))
  end do
  call Grid % Save_Debug_Vtu("a_diagonal", inside_name="a_diagonal", inside_cell=work)
# endif

  call Profiler % Stop('Form_Pressure_Matrix')

  end subroutine
