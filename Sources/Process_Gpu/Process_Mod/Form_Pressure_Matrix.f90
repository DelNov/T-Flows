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
  type(Grid_Type),       pointer :: Grid
  type(Sparse_Con_Type), pointer :: Acon
  type(Sparse_Val_Type), pointer :: Aval, Mval
  integer                        :: s, c1, c2, c
  real                           :: a12
  real, allocatable              :: work(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Form_Pressure_Matrix')

  ! Take some aliases
  Grid => Flow % pnt_grid
  Acon => Flow % Nat % C
  Aval => Flow % Nat % A
  Mval => Flow % Nat % M

  !----------------------------------------------!
  !   This will create race conditions on GPUs   !
  !----------------------------------------------!
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    Assert(c1 .gt. 0)
    Assert(c2 .gt. 0)

    ! Calculate coeficients for the pressure matrix
    ! Units: m * m^3 s / kg = m^4 s / kg
    a12 = Acon % fc(s) * 0.5 * (Mval % v_m(c1) + Mval % v_m(c2))
    Aval % val(Acon % pos(1,s)) = -a12
    Aval % val(Acon % pos(2,s)) = -a12
    Aval % val(Acon % dia(c1))  = Aval % val(Acon % dia(c1)) + a12
    Aval % val(Acon % dia(c2))  = Aval % val(Acon % dia(c2)) + a12

  end do

  ! De-singularize the system matrix ... just like this, ad-hoc
  do c = Cells_In_Domain()
    Aval % val(Acon % dia(c)) = Aval % val(Acon % dia(c)) * (1.0 + MICRO)
  end do

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = 1, Grid % n_cells
    work(c) = Aval % val(Acon % dia(c))
  end do
  call Grid % Save_Debug_Vtu("a_diagonal",              &
                             inside_name="a_diagonal",  &
                             inside_cell=work)
# endif

  call Profiler % Stop('Form_Pressure_Matrix')

  end subroutine
