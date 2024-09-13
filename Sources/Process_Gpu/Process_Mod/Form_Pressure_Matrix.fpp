!==============================================================================!
  subroutine Form_Pressure_Matrix(Process, Acon, Aval, Flow, Grid)
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
  class(Process_Type)                 :: Process
  type(Sparse_Con_Type),       target :: Acon
  type(Sparse_Val_Type),       target :: Aval
  type(Field_Type),            target :: Flow
  type(Grid_Type), intent(in), target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  integer                        :: s, c1, c2, c, i_cel, i, nz
  real                           :: a12
  real, allocatable              :: work(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Pressure_Matrix')

  if(.not. Flow % Nat % use_one_matrix) then
    if(Flow % Nat % A(MATRIX_PP) % formed) return
  end if

  val => Aval % val
  dia => Flow % Nat % C % dia
  pos => Flow % Nat % C % pos
  fc  => Flow % Nat % C % fc
  nz  =  Flow % Nat % C % nonzeros

  !$tf-acc loop begin
  do i = 1, nz
    val(i) = 0.0
  end do
  !$tf-acc loop end

  !--------------------------------------------------------------------!
  !   This is cell based and will not create race conditions on GPUs   !
  !--------------------------------------------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        a12 = fc(s) * Face_Value(s, Flow % v_m(c1), Flow % v_m(c2))
        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a12
        end if
        val(dia(c1)) = val(dia(c1)) + a12
      end if
    end do

  end do
  !$tf-acc loop end

  ! De-singularize the system matrix ... just like this, ad-hoc
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    val(dia(c)) = val(dia(c)) * (1.0 + MICRO)
  end do
  !$tf-acc loop end

  !-------------------------------!
  !   Mark the matrix as formed   !
  !-------------------------------!
  Aval % formed = .true.

# if T_FLOWS_DEBUG == 1
  allocate(work(Grid % n_cells));  work(:) = 0.0
  do c = 1, Grid % n_cells  ! this is for debugging and should be on CPU
    work(c) = val(dia(c))
  end do
  call Grid % Save_Debug_Vtu("a_diagonal",              &
                             inside_name="a_diagonal",  &
                             inside_cell=work)
# endif

  call Profiler % Stop('Form_Pressure_Matrix')

  end subroutine
