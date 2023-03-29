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
  type(Elem_Type),  pointer :: nx, ny, nz
  integer                   :: c, i_fac, fu, e
  real                      :: t_s, x, y, z, r_x, r_y, r_z
  real                      :: radius, dist, dx, dy, dz, xi
  real                      :: beta, t_inf, deltat, thick, ttmp
  real                      :: coef1, coef2, coef3, coef4, coef5, grad
  real, pointer, contiguous :: var(:), res(:), f_int(:)
!  integer                   :: c, c1, c2, s, i_probe, c_inters, n, i_ele, g, l
!  integer                   :: i, j, fu, e, ee, n_cylinders, i_cell, cn
!  real                      :: fs, min_dist, dist, glo_dist, t_s, t_liq
!  real                      :: radius, height, x, y, z, lambda
!  real                      :: p1_x, p1_y, p1_z, distance, cond_1, cond_2
!  real                      :: phi_c1, phi_c2, tot_gradT, tmp 
  !real,         allocatable :: nx(:), ny(:), nz(:)
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  t    => Flow % t
  t_s = Vof % t_sat
  
!   deltat = 1.25
!   t_inf = t_s + deltat
  r_x = 0.0
  r_y = 0.0
  r_z = 0.0

  ! READ THE RADIUS YOU DESIRE
  CALL FILE % OPEN_FOR_READING_ASCII('DESIRED_R', FU)
  READ(FU, *) radius 
  CLOSE(FU)
  PRINT *, 'DESIRED RADIUS:              : ', radius
  
    do c = 1, Grid % n_cells !- Grid % Comm % n_buff_cells
     x = Grid % xc(c)
     y = Grid % yc(c)
     z = Grid % zc(c)
     dist = sqrt((x-r_x)**2 + (y-r_y)**2 + (z-r_z)**2)
     if(fun % n(c) .lt. 0.5) then
       if(Flow % heat_transfer) then
         Flow % t % n(c) = t_s 
         Flow % t % o(c) = t_s 
       end if 
     else ! temperature gradient = 50
       if(Flow % heat_transfer) then
         Flow % t % n(c) = t_s + 10000*(dist-radius)
         Flow % t % o(c) = t_s + 10000*(dist-radius)
       end if
     end if
  end do

  if(DEBUG) then
     
!     call Grid % Save_Debug_Vtu("grad-t-init",                          &
!                                scalar_cell = t % n,                    &
!                                scalar_name = "t",                      &
!                                vector_cell = (/t % x, t % y, t % z/),  &
!                                vector_name = "grad-t-init")
!     grad = 0.0
!     do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
!        e = Vof % Front % elem_in_cell(c)
!        if(e > 0) then
!          grad = sqrt((t % x(c))**2 + (t % y(c))**2 + (t % z(c))**2)
!          write(1000,*) c, grad
!        end if
!        
!     end do
!     
!     Gradients
!     if(.not. Flow % mass_transfer) then
!        call Flow % Grad_Variable(t)
! 
!     if mass transfer, compute gradients with saturation temperature
!     else
!        call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
!     end if
!     
!         call Grid % Save_Debug_Vtu("grad-t-init-front",                          &
!                                scalar_cell = t % n,                    &
!                                scalar_name = "t-front",                      &
!                                vector_cell = (/t % x, t % y, t % z/),  &
!                                vector_name = "grad-t-init-front")
!     grad = 0.0
!     do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
!        e = Vof % Front % elem_in_cell(c)
!        if(e > 0) then
!          grad = sqrt((t % x(c))**2 + (t % y(c))**2 + (t % z(c))**2)
!          write(1001,*) c, grad
!        end if
!        
!     end do
  end if

!   stop
  call Flow % Grad_Variable(t)
  
  end subroutine

