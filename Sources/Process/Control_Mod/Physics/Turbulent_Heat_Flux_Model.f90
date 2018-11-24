!==============================================================================!
  subroutine Control_Mod_Turbulent_Heat_Flux_Model(verbose)
!------------------------------------------------------------------------------!
!   Reading turbulent heat flux model from the control file.                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod,       only: this_proc, Comm_Mod_End
  use Turbulence_Mod, only: turbulent_heat_flux_model,    &
                            SGDH,                         &
                            GGDH,                         &   
                            AFM 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENT_HEAT_FLUX_MODEL', 'SGDH',  &
                                   val, verbose)
  call To_Upper_Case(val)

  !-----------------------------!
  !   Select turbulence model   !
  !-----------------------------!
  select case(val)

    case('SGDH')
      turbulent_heat_flux_model = SGDH
    case('GGDH')
      turbulent_heat_flux_model = GGDH
    case('AFM')
      turbulent_heat_flux_model = AFM

    case default
      if(this_proc < 2) then
        print *, '# ERROR!  Unknown turbulent heat flux model :', trim(val)
        print *, '# Exiting!'
      end if
      call Comm_Mod_End

  end select

  end subroutine
