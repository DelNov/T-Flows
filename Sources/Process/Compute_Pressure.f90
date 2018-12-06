!==============================================================================!
  subroutine Compute_Pressure(flow, sol, dt, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the SIMPLE method.                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Comm_Mod
  use Const_Mod
  use Grid_Mod,    only: Grid_Type
  use Bulk_Mod,    only: Bulk_Type
  use Info_Mod,    only: Info_Mod_Iter_Fill_At
  use Solver_Mod,  only: Solver_Type, Bicg, Cg, Cgs, Acm
  use Matrix_Mod,  only: Matrix_Type
  use Control_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
  real                      :: dt
  integer                   :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  real,              pointer :: flux(:)
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, niter
  real                       :: u_f, v_f, w_f, a12, fs
  real                       :: ini_res, tol, mass_err
  real                       :: px_f, py_f, pz_f
  character(len=80)          :: solver, precond
  real                       :: urf           ! under-relaxation factor
  real                       :: p_max, p_min, p_nor, p_nor_c
!==============================================================================!
!
!   The form of equations which I am solving:
!
!      /               /
!     |               |
!     | rho u dS = dt | GRAD pp dS
!     |               |
!    /               /
!
!   Dimension of the system under consideration
!
!     [App] {pp} = {bpp}               [kg/s]
!
!   Dimensions of certain variables
!
!     app            [ms]
!     pp             [kg/ms^2]
!     b              [kg/s]
!     flux           [kg/s]
!
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  flux => flow % flux
  u    => flow % u
  v    => flow % v
  w    => flow % w
  p    => flow % p
  pp   => flow % pp
  a    => sol % a
  b    => sol % b % val

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(flow, dt, ini)

  !--------------------------------------------------!
  !   Find the value for normalization of pressure   !
  !--------------------------------------------------!

  ! From control file
  call Control_Mod_Normalization_For_Pressure_Solver(p_nor_c)

  ! Calculate pressure magnitude for normalization of pressure solution
  p_max = -HUGE
  p_min = +HUGE
  do c = 1, grid % n_cells
    p_max = max(p_max, p % n(c))
    p_min = min(p_min, p % n(c))
  end do
  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  ! Normalize pressure with the maximum of pressure difference, 
  ! value defined in control file and pressure drops.
  p_nor = max( (p_max-p_min), p_nor_c, abs(bulk % p_drop_x),  &
                                       abs(bulk % p_drop_y),  &
                                       abs(bulk % p_drop_z) )

  ! Initialize matrix and right hand side
  b       = 0.0 
  a % val = 0.0

  !-----------------------------------------!
  !   Initialize the pressure corrections   !
  !-----------------------------------------!
  pp % n = 0.0 

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if(c2 > 0) then

      ! Interpolate velocity 
      u_f = fs * u % n(c1) + (1.0-fs) * u % n(c2)
      v_f = fs * v % n(c1) + (1.0-fs) * v % n(c2)
      w_f = fs * w % n(c1) + (1.0-fs) * w % n(c2)

      ! Calculate coeficients for the system matrix
      if(c2 > 0) then
        a12 = 0.5 * density * a % fc(s) *        &
           (  grid % vol(c1) / a % sav(c1)       &
            + grid % vol(c2) / a % sav(c2) )
        a % val(a % pos(1,s)) = -a12
        a % val(a % pos(2,s)) = -a12
        a % val(a % dia(c1))  = a % val(a % dia(c1)) +  a12
        a % val(a % dia(c2))  = a % val(a % dia(c2)) +  a12
      else  ! I am somewhat surprised this part is here
        a12 = 0.5 * density * a % fc(s) *        &
             (  grid % vol(c1) / a % sav(c1)     &
              + grid % vol(c2) / a % sav(c2) )
        a % val(a % pos(1,s)) = -a12
        a % val(a % dia(c1))  = a % val(a % dia(c1)) +  a12
      end if 

      ! Interpolate pressure gradients
      px_f = 0.5*( p % x(c1) + p % x(c2) ) * grid % dx(s)
      py_f = 0.5*( p % y(c1) + p % y(c2) ) * grid % dy(s)
      pz_f = 0.5*( p % z(c1) + p % z(c2) ) * grid % dz(s)

      ! Calculate flux through cell face
      flux(s) = density * (  u_f*grid % sx(s)       &
                           + v_f*grid % sy(s)       &
                           + w_f*grid % sz(s) )     &
              + a12 * (p % n(c1) - p % n(c2))       &
              + a12 * (px_f + py_f + pz_f)

      b(c1)=b(c1)-flux(s)
      if(c2  > 0) b(c2)=b(c2)+flux(s)

    ! Side is on the boundary
    else ! (c2 < 0)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        u_f = u % n(c2)
        v_f = v % n(c2)
        w_f = w % n(c2)
        flux(s) = density * (  u_f * grid % sx(s)  &
                             + v_f * grid % sy(s)  &
                             + w_f * grid % sz(s) )
        b(c1) = b(c1)-flux(s)
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW .or.   &
              Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) then
        u_f = u % n(c2)
        v_f = v % n(c2)
        w_f = w % n(c2)
        flux(s) = density * (  u_f*grid % sx(s)  &
                             + v_f*grid % sy(s)  &
                             + w_f*grid % sz(s) )
        b(c1) = b(c1)-flux(s)

        a12 = density * a % fc(s) * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then
        u_f = u % n(c1)
        v_f = v % n(c1)
        w_f = w % n(c1)
        flux(s) = density * (  u_f * grid % sx(s)  &
                             + v_f * grid % sy(s)  &
                             + w_f * grid % sz(s) )
        b(c1) = b(c1)-flux(s)

        a12 = density * a % fc(s) * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12

      else  ! it is SYMMETRY
        flux(s) = 0.0
      end if
    end if

  end do

  mass_err = 0.0
  do c = 1, grid % n_cells
    mass_err = max(mass_err, abs(b(c)))
  end do

  ! Don't solve the pressure corection too accurate.
  ! Value 1.e-18 blows the solution.
  ! Value 1.e-12 keeps the solution stable
  call Control_Mod_Tolerance_For_Pressure_Solver(tol)

  ! Get solver and matrix precondioner
  call Control_Mod_Solver_For_Pressure(solver)
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set the default value for number of iterations / v-cycles
  niter = 40
  call Control_Mod_Max_Iterations_For_Pressure_Solver(niter)

  if(solver .eq. 'ACM') then
    tol   = PICO
    call Acm(sol,       &
             pp % n,    &
             b,         &
             precond,   &
             niter,     &     ! number of V cycles
             tol,       &
             ini_res,   &
             pp % res,  &
             norm = p_nor)    ! last argument: number for normalisation
    stop
  else
    call Cg(sol,       &
            pp % n,    &
            b,         &
            precond,   &
            niter,     &      ! number of iterations
            tol,       &
            ini_res,   &
            pp % res,  &
            norm = p_nor)     ! last argument: number for normalisation
  end if

  call Info_Mod_Iter_Fill_At(1, 3, pp % name, niter, pp % res)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  call Control_Mod_Simple_Underrelaxation_For_Pressure(urf)  ! retreive urf
  do c = 1, grid % n_cells
    p % n(c)  =  p % n(c)  +  urf  *  pp % n(c)
  end do

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  p_max  = maxval(p % n(1:grid % n_cells))
  p_min  = minval(p % n(1:grid % n_cells))

  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)

  p % n(:) = p % n(:) - 0.5*(p_max+p_min)

  call Comm_Mod_Exchange_Real(grid, pp % n)

  ! User function
  call User_Mod_End_Of_Compute_Pressure(flow, dt, ini)

  end subroutine
