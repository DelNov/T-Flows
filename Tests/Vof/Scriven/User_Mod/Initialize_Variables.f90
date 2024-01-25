include './User_Mod/Vof_Initialization_Ellipsoid.f90'

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun, t
  type(Vof_Type),   pointer :: t_sat, latent_heat
  integer                   :: c, fu
  real                      :: t_s, x, y, z, r_x, r_y, r_z
  real                      :: dist
  real                      :: t_inf
  real, allocatable         :: r_in(:), tpr_in(:)
  real                      :: r_min, r_max, tmp1, tmp2
  integer                   :: i, ndata
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  t    => Flow % t
  t_s = Vof % t_sat

  ! VOF
  ! Initialize the whole domain as 0.0
  do c = 1, Grid % n_cells
    fun % n(c) = 0.0
  end do

  ! Ellipsoid:
  !WRITE(*,*)"call Vof_Initialization_Ellipsoid"
  call Vof_Initialization_Ellipsoid(Vof)

  call Grid % Exchange_Cells_Real(fun % n)

  ! Old value
  fun % o(:) = fun % n(:)

  ! Temperature
  call File % Open_For_Reading_Ascii("dT1.25_R50microns.txt", fu)
    READ(fu,*)      !header
    DO i = 1,10000
      READ(fu,*,end=999,err=999)tmp1,tmp2
    ENDDO
 999 CONTINUE
    ndata = i-1
    ALLOCATE(r_in(ndata),tpr_in(ndata))

    rewind(fu)
    READ(fu,*)     !header
    DO i = 1, ndata
      READ(fu,*)r_in(i),tpr_in(i)
      !WRITE(*,*)r_in(i),tpr_in(i)
    ENDDO
  close(fu)

  t_inf = tpr_in(ndata)
  r_x = 0.0
  r_y = 0.0
  r_z = 0.0
  r_min = r_in(1)
  r_max = r_in(ndata)

  ! Initialize temperatures from file
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    x = Grid % xc(c)
    y = Grid % yc(c)
    z = Grid % zc(c)
    dist = sqrt((x-r_x)**2 + (y-r_y)**2 + (z-r_z)**2)
    if (dist < r_min) then
      Flow % t % n(c) = t_s 
      Flow % t % o(c) = t_s
    else if (dist > r_max) then
      Flow % t % n(c) = t_inf
      Flow % t % o(c) = t_inf
    else
      Flow % t % n(c) = tpr_interpolate(dist,r_in,tpr_in,ndata)
      Flow % t % o(c) = tpr_interpolate(dist,r_in,tpr_in,ndata)
    end if
  end do

  DEALLOCATE(r_in, tpr_in)
  
  call Flow % Grad_Variable(t)
  
  end subroutine

  !SUBROUTINE tpr_interpolate(dist,r_in,tpr_in,ndata)
  REAL FUNCTION tpr_interpolate(dist,r_in,tpr_in,ndata)
    INTEGER::ndata,i
    REAL::dist,r_in(ndata),tpr_in(ndata),ratio
    DO i = 1, ndata-1
      if (r_in(i)<=dist .AND. dist<r_in(i+1)<dist) then
        ratio = (dist-r_in(i))/(r_in(i+1)-r_in(i))
        tpr_interpolate = tpr_in(i)*(1.0-ratio) + tpr_in(i+1)*ratio
        CYCLE
      endif
    ENDDO
    !WRITE(*,*)dist,tpr_interpolate
  END FUNCTION tpr_interpolate
