#include "./Vof_Initialization_Ellipsoid.f90"

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun
  real,             pointer :: dt
  integer                   :: c,s,c2
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  dt   => Flow % dt

  ! Initialize the whole domain as 0.0
  do c = Cells_In_Domain_And_Buffers()
    fun % n(c) = 0.0
  end do

  ! Ellipsoid:
  write(*,*) "call Vof_Initialization_Ellipsoid"
  call Vof_Initialization_Ellipsoid(Vof)

  call Grid % Exchange_Cells_Real(fun % n)

  ! Old value
  fun % o(:) = fun % n(:)

  !----------------------------------------!
  !   Correct for contact angle at walls   !
  !----------------------------------------!
  !do s = 1, Grid % n_faces
  ! c2 = Grid % faces_c(2,s)
  !  if(c2 < 0) then
  !    if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.    &
  !     Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then
  !      WRITE(*,*)"Curvature",fun % q(c2)
  !      fun % q(c2) = 180.0
  !      !fun % q(c2) = 90.0
  !    endif
  !  endif
  !enddo

  end subroutine
