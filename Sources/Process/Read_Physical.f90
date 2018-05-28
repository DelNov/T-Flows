!==============================================================================!
  subroutine Read_Physical(grid, restar)
!------------------------------------------------------------------------------!
!   Reads details about physial models.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Rans_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: restar
!-----------------------------------[Locals]-----------------------------------!
  integer :: m
!==============================================================================!

  if(this_proc < 2) &
    print *, 'Calling for gravitationa vector'
  call Control_Mod_Gravitational_Vector(grav_x, grav_y, grav_z)

  if(this_proc < 2) &
    print *, 'Calling for turbulence model'
  call Control_Mod_Turbulence_Model(.true.)

  call Control_Mod_Rough_Walls(.true.)

  if(turbulence_model == K_EPS) then
    call Constants_K_Eps()
  endif

  if(turbulence_model == REYNOLDS_STRESS) then
    call Constants_Reynolds_Stress()
  endif

  if(turbulence_model == HANJALIC_JAKIRLIC) then
    call Constants_Hanjalic_Jakirlic()
  endif

  if(turbulence_model == K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    call Constants_K_Eps_Zeta_F()
  endif

  if(turbulence_model == SPALART_ALLMARAS .or.  &
     turbulence_model == DES_SPALART) then
    call Constants_Spalart_Allmaras()
  end if

  ! Wall velocity
  if(.not. restar) then
    do m=1,grid % n_materials
      call Control_Mod_Mass_Flow_Rates(bulk(m) % p_drop_x,  &
                                       bulk(m) % p_drop_y,  &
                                       bulk(m) % p_drop_z)
    end do
  end if

  ! Mass fluxes
  if(.not. restar) then
    do m=1,grid % n_materials
      call Control_Mod_Mass_Flow_Rates(bulk(m) % flux_x_o,  &
                                       bulk(m) % flux_y_o,  &
                                       bulk(m) % flux_z_o)
    end do
  end if

  end subroutine
