!==============================================================================!
  subroutine Physical_Properties(Rc, Grid, Flow)
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
  type(Field_Type)                      :: Flow  !! flow object
!-----------------------------------[Locals]-----------------------------------!
  real :: dens_const, visc_const
  real :: capa_const, cond_const
  real :: diff_const
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Give some sign
  O_Print '(a)', ' # Reading about physical properties'

  ! Read constant (defualt) values
  call Control % Dynamic_Viscosity(visc_const)
  call Control % Mass_Density     (dens_const)
  if(Flow % heat_transfer) then
    call Control % Heat_Capacity       (capa_const)
    call Control % Thermal_Conductivity(cond_const)
  end if
  if (Flow % n_scalars .gt. 0) then
    call Control % Scalars_Diffusivity (diff_const)
  end if

  Flow % density  (:) = dens_const
  Flow % viscosity(:) = visc_const
  if(Flow % heat_transfer) then
    Flow % capacity    (:) = capa_const
    Flow % conductivity(:) = cond_const
  end if
  if(Flow % n_scalars .gt. 0) then
    Flow % diffusivity (:) = diff_const
  end if

  end subroutine
