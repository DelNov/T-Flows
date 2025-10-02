!==============================================================================!
  subroutine Numerical_Schemes(Rc, Grid, Flow, Turb)
!------------------------------------------------------------------------------!
!>  This is s a simplified version from the same subroutine in Process_Cpu
!>  as it reads only boundary conditions releated to momentum and enthalpy
!>  conservation equations.  Hopefully, as more modules are ported to
!>  Process_Gpu, this source file will get closer and closer to its origin
!>  from Process_Cpu.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc    !! parent class
  type(Grid_Type)                       :: Grid  !! grid object
  type(Field_Type), target              :: Flow  !! flow object
  type(Turb_Type),  target              :: Turb  !! turbulence object
!----------------------------------[Locals]------------------------------------!
  type(Var_Type), pointer :: tq, ui, phi
  character(SL)           :: name
  integer                 :: i, sc
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Give some sign
  O_Print '(a)', ' # Reading info about discretization schemes'

  !------------------------------------------!
  !   Pressure velocity coupling algorithm   !
  !------------------------------------------!

  ! Report volume balance (in a separate file)
  call Control % Report_Volume_Balance(Flow % rep_vol_balance, .false.)

  !----------------------------------!
  !   Gradient computation methods   !
  !----------------------------------!

  ! Tolerance and max iterations for computation of gradients
  call Control % Max_Least_Squares_Gradients_Iterations(Flow % least_miter,  &
                                                          .false.)
  !-------------------------!
  !   Related to momentum   !
  !-------------------------!
  do i = 1, 3
    if(i .eq. 1) ui => Flow % u
    if(i .eq. 2) ui => Flow % v
    if(i .eq. 3) ui => Flow % w
    call Control % Blending_Coefficient_For_Momentum  (ui % blend)
    call Control % Simple_Underrelaxation_For_Momentum(ui % urf)
    call Control % Blend_System_Matrices(ui % blend_matrix, .false.)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  call Control % Simple_Underrelaxation_For_Pressure(Flow % pp % urf)
  call Control % Blend_System_Matrices(Flow % pp % blend_matrix, .false.)

  !------------------------------!
  !   Related to wall distance   ! (nothing yet)
  !------------------------------!

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(Flow % heat_transfer) then
    call Control % Blending_Coefficient_For_Energy  (Flow % t % blend)
    call Control % Simple_Underrelaxation_For_Energy(Flow % t % urf)
    call Control % Blend_System_Matrices(Flow % t % blend_matrix, .false.)
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Control % Blending_Coefficient_For_Scalars  (phi % blend)
    call Control % Simple_Underrelaxation_For_Scalars(phi % urf)
    call Control % Blend_System_Matrices(phi % blend_matrix, .false.)
  end do

  !------------------------------!
  !   All turbulent quantities   !
  !------------------------------!
  do i = 1, 12
    nullify(tq)
    if(i .eq. 5) tq => Turb % vis
    if(associated(tq)) then
      call Control % Blending_Coefficient_For_Turbulence  (tq % blend)
      call Control % Simple_Underrelaxation_For_Turbulence(tq % urf)
      call Control % Blend_System_Matrices(tq % blend_matrix, .false.)
    end if
  end do

  end subroutine
