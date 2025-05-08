!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Grid, Flow, Turb)
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: Grid
  type(Field_Type) :: Flow
  type(Turb_Type)  :: Turb
!------------------------------[Local parameters]------------------------------!
  real, parameter :: R_MAX = 0.0025  ! m
  real, parameter :: U_MAX = 1.0   ! m/s
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2, reg
  real    :: r, u
!------------------------[Avoid unused parent warning]-------------------------!
!==============================================================================!

  !--------------------------------------------------------!
  !                                                        !
  !   Initialize velocities and volume fluxes everywhere   !
  !                                                        !
  !--------------------------------------------------------!

  !----------------------------------------------------------!
  !   Initialize velocities in all cells inside the domain   !
  !----------------------------------------------------------!
  O_Print '(a)', ' # Prescribing parabolic velocity profile'
  do c = Cells_In_Domain()
    ! Pipe section is parallel to yz plane
    r = sqrt(Grid % yc(c)**2 + Grid % zc(c)**2)
    Flow % u % n(c) = U_MAX * (1 - r**2 / R_MAX**2)
  end do

  !----------------------------------------------!
  !   Initialize volume fluxes at inside faces   !
  !----------------------------------------------!
  !print *, 'Flux before:', Flow % v_flux % n(10945)
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    u = 0.5 * (Flow % u % n(c1) + Flow % u % n(c2))
    Flow % v_flux % n(s) = Grid % sx(s) * u
  end do
  !print *, 'Flux after:', Flow % v_flux % n(10945)

  !-----------------------------------------------------------!
  !   Initialize velocities in all boundary cells and faces   !
  !-----------------------------------------------------------!
  do reg = 1, Grid % n_regions
    if(Grid % region % name(reg) .eq. 'INLET') then
      print *, 'Flux before:', Flow % v_flux % n(9396)
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1, s)
        c2 = Grid % faces_c(2, s)  ! boundary cell
        r = sqrt(Grid % yc(c2)**2 + Grid % zc(c2)**2)
        Flow % u % n(c2) = U_MAX * (1 - r**2 / R_MAX**2)
        Flow % v_flux % n(s) = + Grid % sx(s) * Flow % u % n(c2)
      end do
      print *, 'Flux after:', Flow % v_flux % n(9396)
    end if
    if(Grid % region % name(reg) .eq. 'OUTLET') then
      print *, 'Flux before:', Flow % v_flux % n(10944)
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1, s)
        c2 = Grid % faces_c(2, s)  ! boundary cell
        r = sqrt(Grid % yc(c2)**2 + Grid % zc(c2)**2)
        Flow % u % n(c2) = U_MAX * (1 - r**2 / R_MAX**2)
        Flow % v_flux % n(s) = + Grid % sx(s) * Flow % u % n(c2)
      end do
      print *, 'Flux after:', Flow % v_flux % n(10944)
    end if
  end do

  !-----------------------------------------!
  !                                         !
  !   Allocate memory for reaction rates,   !
  !   initialize values, and copy to GPU    !
  !                                         !
  !-----------------------------------------!
  allocate(rate(6,6))
  rate(:,:) = 0.0

  ! Here we will put some magic from GEMS (if it exists)
  rate(PB_G,  PBI_G) = 0.0e-3
  rate(PB_G,  I_G)   = 0.0e-3
  rate(PB_S,  PB_G)  = 0.0e-3

  call Gpu % Matrix_Real_Copy_To_Device(rate)

  end subroutine
