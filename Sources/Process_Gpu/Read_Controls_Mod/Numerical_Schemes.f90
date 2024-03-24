!==============================================================================!
  subroutine Numerical_Schemes(Rc, Flow)
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
  type(Field_Type), target              :: Flow  !! flow object
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: ui, phi
  character(SL)            :: name
  integer                  :: i, sc
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Give some sign
  O_Print '(a)', ' # Reading info about discretization schemes'

  ! Take alias
  Grid => Flow % pnt_grid

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
    call Control % Blending_Coefficient_For_Momentum    (ui % blend)
    call Control % Simple_Underrelaxation_For_Momentum  (ui % urf)
  end do

  !-------------------------!
  !   Related to pressure   !
  !-------------------------!
  call Control % Simple_Underrelaxation_For_Pressure(Flow % pp % urf)

  !------------------------------!
  !   Related to heat transfer   !
  !------------------------------!
  if(Flow % heat_transfer) then
    call Control % Blending_Coefficient_For_Energy  (Flow % t % blend)
    call Control % Simple_Underrelaxation_For_Energy(Flow % t % urf)
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    call Control % Blending_Coefficient_For_Scalars            (phi % blend)
    call Control % Simple_Underrelaxation_For_Scalars          (phi % urf)
  end do

  end subroutine