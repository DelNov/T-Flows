!==============================================================================!
  subroutine Vis_T_Subgrid(Turb, Flow, Grid)
!------------------------------------------------------------------------------!
!>  Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.           !
!   In case that, in parallel executions, the subdomain does not have          !
!   any nearwall cells, the nearest_wall_cell(c) is zero.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb  !! parent class
  type(Field_Type), target :: Flow
  type(Grid_Type),  target :: Grid
!------------------------------[Local parameters]------------------------------!
  real, parameter :: A_POW = 8.3
  real, parameter :: B_POW = 1.0/7.0
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
  real    :: u_tau, u_tan, cs, nu, lf
!==============================================================================!

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then

    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()
      lf = Grid % vol(c)**ONE_THIRD

      nu = Flow % viscosity(c) / Flow % density(c)

      ! Tangential velocity.  Here it assumed that, as you approach the
      ! wall, the tangential velocity component is dominant and that the
      ! magnitude of velocity is close to tangential component.
      u_tan = sqrt(  Flow % u % n(c)**2   &
                   + Flow % v % n(c)**2   &
                   + Flow % w % n(c)**2)

      ! Calculate u_tau, y+ and perform Van Driest damping
      u_tau = (u_tan/A_POW * (nu/Grid % wall_dist(c))**B_POW)   &
              ** (1.0/(1.0+B_POW))
      Turb % y_plus(c) = Grid % wall_dist(c) * u_tau / Flow % viscosity(c)
      cs = Turb % c_smag * (1.0 - exp(-Turb % y_plus(c) / 25.0))

      Turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (cs*cs)            &  ! cs^2
                      * Flow % shear(c)
    end do
    !$tf-acc loop end

  else if(Turb % model .eq. LES_WALE) then

    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()
      lf = Grid % vol(c)**ONE_THIRD
      Turb % vis_t(c) = Flow % density(c)  &
                      * (lf*lf)            &  ! delta^2
                      * (0.5*0.5)          &  ! cs^2
                      * Turb % wale_v(c)
    end do
    !$tf-acc loop end

  end if

  !-------------------!
  !   Wall function   !
  !-------------------!
  call Turb % Wall_Function(Grid, Flow)

  ! Refresh buffers
  call Grid % Exchange_Cells_Real(Turb % vis_t)

  end subroutine
