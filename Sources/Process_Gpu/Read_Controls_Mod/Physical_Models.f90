!==============================================================================!
  subroutine Physical_Models(Rc, Flow)
!------------------------------------------------------------------------------!
!>  This is s a simplified version from the same subroutine in Process_Cpu
!>  as it reads only boundary conditions releated to momentum and enthalpy
!>  conservation equations.  Hopefully, as more modules are ported to
!>  Process_Gpu, this source file will get closer and closer to its origin
!>  from Process_Cpu.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc     !! parent class
  type(Field_Type), target              :: Flow   !! flow object
!----------------------------------[Locals]------------------------------------!
  type(Bulk_Type), pointer :: bulk
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Give some sign
  if(First_Proc())  &
    print '(a)', ' # Reading about physical models'

  ! Take aliases
  bulk => Flow % bulk

  !-------------------------------------------!
  !                                           !
  !   Related to heat transfer and bouyancy   !
  !                                           !
  !-------------------------------------------!
  call Control % Heat_Transfer(Flow % heat_transfer, verbose = .true.)

  !------------------------------------!
  !                                    !
  !   Pressure drops and mass fluxes   !
  !                                    !
  !------------------------------------!
  call Control % Pressure_Drops(bulk % p_drop_x,  &
                                bulk % p_drop_y,  &
                                bulk % p_drop_z)
  call Control % Bulk_Velocities(bulk % u_o,  &
                                 bulk % v_o,  &
                                 bulk % w_o)

  !-----------------------!
  !                       !
  !   Number of scalars   !
  !                       !
  !-----------------------!
  call Control % Number_Of_Scalars(Flow % n_scalars, verbose = .true.)

  end subroutine
