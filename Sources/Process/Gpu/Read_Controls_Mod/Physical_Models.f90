!==============================================================================!
  subroutine Physical_Models(Rc, Flow, Turb)
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
  type(Turb_Type),  target              :: Turb   !! turbulence object
!----------------------------------[Locals]------------------------------------!
  type(Bulk_Type), pointer :: bulk
  character(SL)            :: name
  logical                  :: found
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
  call Control % Gravitational_Vector(Flow % grav_x,  &
                                      Flow % grav_y,  &
                                      Flow % grav_z, .true.)
  call Control % Buoyancy(name, .true.)
  select case(name)
    case('NONE')
      Flow % buoyancy = NO_BUOYANCY
    case('DENSITY')
      Flow % buoyancy = DENSITY_DRIVEN
    case('THERMAL')
      Flow % buoyancy = THERMALLY_DRIVEN
    case default
      call Message % Error(60,                                       &
                           'Unknown buoyancy model: '//trim(name)//  &
                           '.  Exiting!',                            &
                           file=__FILE__, line=__LINE__)
  end select

  call Control % Reference_Density           (Flow % dens_ref, .true.)
  call Control % Reference_Temperature       (Flow % t_ref,    .true.)
  call Control % Volume_Expansion_Coefficient(Flow % beta,     .true.)

  !---------------------------!
  !                           !
  !   Related to turbulence   !
  !                           !
  !---------------------------!
  call Control % Turbulence_Model(name, .true.)
  select case(name)

    case('NONE')
      Turb % model = NO_TURBULENCE_MODEL
    case('LES_SMAGORINSKY')
      Turb % model = LES_SMAGORINSKY
    case('LES_WALE')
      Turb % model = LES_WALE
    case('SPALART_ALLMARAS')
      Turb % model = SPALART_ALLMARAS

    case default
      call Message % Error(60,                                         &
                           'Unknown turbulence model: '//trim(name)//  &
                           '.  Exiting!',                              &
                           file=__FILE__, line=__LINE__)
  end select

  !-------------------------------------------------------------------------!
  !   Initialization of model constants depending on the turbulence model   !
  !-------------------------------------------------------------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Control % Smagorinsky_Constant(Turb % c_smag, .true.)
  end if


  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Turb % Const_Spalart_Allmaras()
  end if

  if(Turb % model .eq. LES_DYNAMIC     .or.  &
     Turb % model .eq. LES_SMAGORINSKY .or.  &
     Turb % model .eq. LES_WALE) then
    call Turb % Const_Les()
  end if

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
  call Control % Position_At_One_Key('VOLUME_FLOW_RATES', found)
  if(found) then
    call Message % Error(88,                                                  &
      'You specified volume flow rates using the keyword '                 // &
      '"VOLUME_FLOW_RATES" in the control file, which was common '         // &
      'practice in the CPU version of the code.  However, this feature '   // &
      'has been superseded by bulk velocities in the GPU version of the '  // &
      'code, and should be defined using the keyword "BULK_VELOCITIES". '  // &
      ' \n \n Please update the control file accordingly and try running ' // &
      'the simulation again.', file=__FILE__, line=__LINE__)
  end if

  !-----------------------!
  !                       !
  !   Number of scalars   !
  !                       !
  !-----------------------!
  call Control % Number_Of_Scalars(Flow % n_scalars, verbose = .true.)

  end subroutine
