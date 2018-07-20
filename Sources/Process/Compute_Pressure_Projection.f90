!==============================================================================!
  subroutine Compute_Pressure_Fractional(grid, dt, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the fractional step method.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Grid_Mod,     only: Grid_Type  
  use Comm_Mod
  use Info_Mod
  use Solvers_Mod,  only: Bicg, Cg, Cgs
  use Control_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: dt
  integer         :: ini
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter
  real              :: p_max, p_min
  real              :: ini_res, tol, mass_err
  real              :: u_f, v_f, w_f, fs
  real              :: a12
  character(len=80) :: precond
  real              :: urf           ! under-relaxation factor                 
!==============================================================================!
!     
!   The form of equations which are being solved:
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
!     pp,            [kg/ms^2]
!     b              [kg/s]
!     flux           [kg/s]
!   
!------------------------------------------------------------------------------!

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(grid, dt, ini)

  ! Initialize matrix and source term
  a % val = 0.0
  b = 0.0 

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if( c2 > 0 .or.  &
        c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then 

      ! Extract the "centred" pressure terms from cell velocities
      u_f = fs*      (u % n(c1) + p % x(c1)*grid % vol(c1)/a % sav(c1))       &
          + (1.0-fs)*(u % n(c2) + p % x(c2)*grid % vol(c2)/a % sav(c2))

      v_f = fs*      (v % n(c1) + p % y(c1)*grid % vol(c1)/a % sav(c1))       &
          + (1.0-fs)*(v % n(c2) + p % y(c2)*grid % vol(c2)/a % sav(c2))

      w_f = fs*      (w % n(c1) + p % z(c1)*grid % vol(c1)/a % sav(c1))       &
          + (1.0-fs)*(w % n(c2) + p % z(c2)*grid % vol(c2)/a % sav(c2))

      ! Add the "staggered" pressure terms to face velocities
      u_f = u_f + (p % n(c1) - p % n(c2))  &
                * grid % sx(s) * ( fs / a % sav(c1) + (1.0-fs) / a % sav(c2) )
      v_f = v_f + (p % n(c1) - p % n(c2))  &
                * grid % sy(s) * ( fs / a % sav(c1) + (1.0-fs) / a % sav(c2) )
      w_f = w_f + (p % n(c1) - p % n(c2))  &
                * grid % sz(s) * ( fs / a % sav(c1) + (1.0-fs) / a % sav(c2) )

      ! Now calculate the flux through cell face
      flux(s) = density * ( u_f * grid % sx(s) +  &
                            v_f * grid % sy(s) +  &
                            w_f * grid % sz(s) )

      a12 = density * (  grid % sx(s) * grid % sx(s)  &
                       + grid % sy(s) * grid % sy(s)  &
                       + grid % sz(s) * grid % sz(s))
      a12 = a12 * (fs/a % sav(c1) + (1.-fs)/a % sav(c2))

      if(c2  > 0) then 
        a % val(a % pos(1,s)) = -a12
        a % val(a % pos(2,s)) = -a12
        a % val(a % dia(c1))  = a % val(a % dia(c1)) +  a12
        a % val(a % dia(c2))  = a % val(a % dia(c2)) +  a12
      else
        a % bou(c2) = -a12
        a % val(a % dia(c1)) = a % val(a % dia(c1)) +  a12
      endif

      b(c1)=b(c1)-flux(s)
      if(c2  > 0) b(c2)=b(c2)+flux(s)

    ! Face is on the boundary
    else
      u_f = u % n(c2)
      v_f = v % n(c2)
      w_f = w % n(c2)

      flux(s) = density * (  u_f * grid % sx(s)  &
                           + v_f * grid % sy(s)  &
                           + w_f * grid % sz(s) )

      b(c1) = b(c1)-flux(s)
    end if

  end do

  !----------------------------------------!
  !   Initialize the pressure correction   !
  !----------------------------------------!
  pp % n = 0.0 

  mass_err = 0.0
  do c = 1, grid % n_cells
    mass_err = mass_err + abs(b(c))
  end do
  call Comm_Mod_Global_Sum_Real(mass_err)                       

  !--------------------------------------------!
  !   Solve the pressure correction equation   !
  !--------------------------------------------!

  ! Don't solve the pressure corection too accurate.
  ! Value 1.e-18 blows the solution.
  ! Value 1.e-12 keeps the solution stable
  call Control_Mod_Tolerance_For_Pressure_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set the default value for number of iterations
  niter = 200

  ! Over-ride if specified in control file
  call Control_Mod_Max_Iterations_For_Pressure_Solver(niter)

  call Cg(a, pp % n, b, precond, niter, tol, ini_res, pp % res)

  call Info_Mod_Iter_Fill_At(1, 3, pp % name, niter, pp % res)   

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  urf = 1.0

  p % n(:)  =  p % n(:)   +  urf  *  pp % n(:) 

  !----------------------------------!
  !   Normalize the pressure field   !
  !----------------------------------!
  p_max  = maxval(p % n(1:grid % n_cells))
  p_min  = minval(p % n(1:grid % n_cells))

  call Comm_Mod_Global_Max_Real(p_max) 
  call Comm_Mod_Global_Min_Real(p_min) 

  p % n(:)   =  p % n(:)   -  0.5 * (p_max + p_min)

  call Comm_Mod_Exchange(grid, pp % n) 

  ! User function
  call User_Mod_End_Of_Compute_Pressure(grid, dt, ini)

  end subroutine