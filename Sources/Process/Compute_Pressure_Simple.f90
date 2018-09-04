!==============================================================================!
  subroutine Compute_Pressure_Simple(grid, dt, ini)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the S.I.M.p.L.E. method.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Comm_Mod
  use Const_Mod
  use Grid_Mod,     only: Grid_Type
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
  real              :: u_f, v_f, w_f, a12, fs
  real              :: ini_res, tol, mass_err
  real              :: smdpn
  real              :: px_f, py_f, pz_f
  character(len=80) :: precond
  real              :: urf           ! under-relaxation factor                 
  real              :: min_b, max_b
  real              :: p_max, p_min, p_nor
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

  ! User function
  call User_Mod_Beginning_Of_Compute_Pressure(grid, dt, ini)

  ! Calculate velocity magnitude
  p_max = -HUGE
  p_min = +HUGE
  do c = 1, grid % n_cells
    p_max = max(p_max, p % n(c))
    p_min = min(p_min, p % n(c))
  end do
  call Comm_Mod_Global_Max_Real(p_max)
  call Comm_Mod_Global_Min_Real(p_min)
  p_nor = max( (p_max-p_min), MICRO, abs(bulk % p_drop_x),  &
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
    if(c2  > 0 .or.  &
       c2  < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then

      smdpn = (  grid % sx(s)*grid % sx(s)   &
               + grid % sy(s)*grid % sy(s)   &
               + grid % sz(s)*grid % sz(s) ) &
            / (  grid % sx(s)*grid % dx(s)   & 
               + grid % sy(s)*grid % dy(s)   &
               + grid % sz(s)*grid % dz(s) )  

      ! Interpolate velocity 
      u_f = fs * u % n(c1) + (1.0-fs) * u % n(c2)
      v_f = fs * v % n(c1) + (1.0-fs) * v % n(c2)
      w_f = fs * w % n(c1) + (1.0-fs) * w % n(c2)

      ! Calculate coeficients for the system matrix
      if(c2  > 0) then 
        a12 = 0.5 * density * smdpn *            &
           (  grid % vol(c1) / a % sav(c1)       &
            + grid % vol(c2) / a % sav(c2) )
        a % val(a % pos(1,s))  = -a12
        a % val(a % pos(2,s))  = -a12
        a % val(a % dia(c1)) = a % val(a % dia(c1)) +  a12
        a % val(a % dia(c2)) = a % val(a % dia(c2)) +  a12
      else 
        a12 = 0.5 * density * smdpn *            &
             (  grid % vol(c1) / a % sav(c1)     &
              + grid % vol(c2) / a % sav(c2) )
        a % bou(c2)  = -a12
        a % val(a % dia(c1)) = a % val(a % dia(c1)) +  a12
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
        smdpn = (  grid % sx(s) * grid % sx(s)   &
                 + grid % sy(s) * grid % sy(s)   &
                 + grid % sz(s) * grid % sz(s) ) &
              / (  grid % sx(s) * grid % dx(s)   &
                 + grid % sy(s) * grid % dy(s)   &
                 + grid % sz(s) * grid % dz(s) )  
        a12 = density * smdpn * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) +  a12
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) then
        u_f = u % n(c1)
        v_f = v % n(c1)
        w_f = w % n(c1)
        flux(s) = density * (  u_f * grid % sx(s)  &
                             + v_f * grid % sy(s)  &
                             + w_f * grid % sz(s) )
        b(c1) = b(c1)-flux(s)
        smdpn = ( grid % sx(s) * grid % sx(s)   &
                + grid % sy(s) * grid % sy(s)   &
                + grid % sz(s) * grid % sz(s) ) &
              / ( grid % sx(s) * grid % dx(s)   &
                + grid % sy(s) * grid % dy(s)   &
                + grid % sz(s) * grid % dz(s) )
        a12 = density * smdpn * grid % vol(c1) / a % sav(c1)
        a % val(a % dia(c1)) = a % val(a % dia(c1)) +  a12
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

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set the default value for number of iterations
  niter = 40

  ! Over-ride if specified in control file
  call Control_Mod_Max_Iterations_For_Pressure_Solver(niter)

  call Bicg(a,         &
            pp % n,    &
            b,         &
            precond,   &
            niter,     &
            tol,       &
            ini_res,   &
            pp % res,  &
            norm = p_nor)    ! last argument: number for normalisation

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
  ! p_max  = maxval(p % n(1:grid % n_cells))
  ! p_min  = minval(p % n(1:grid % n_cells))

  ! call Comm_Mod_Global_Max_Real(p_max)
  ! call Comm_Mod_Global_Min_Real(p_min)

  ! p % n = p % n - 0.5*(p_max+p_min)

  call Comm_Mod_Exchange_Real(grid, pp % n)

  ! User function
  call User_Mod_End_Of_Compute_Pressure(grid, dt, ini)

  end subroutine
