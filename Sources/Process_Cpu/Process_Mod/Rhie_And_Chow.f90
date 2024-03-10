!==============================================================================!
  subroutine Rhie_And_Chow(Process, Flow, Vof, Nat)
!------------------------------------------------------------------------------!
!>  The Rhie_And_Chow subroutine in T-Flows is employed to calculate face
!>  velocities using the Rhie and Chow interpolation method. This method is
!>  integral for ensuring a non-oscillatory, stable solution for the face
!>  velocities in the computational fluid dynamics simulations.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization and setup: Initiates by establishing necessary variables  !
!     and pointers for the grid, flow field, and other parameters. It also     !
!     sets up the profiler for tracking performance and efficiency.            !
!   * Pressure and surface tension handling: The subroutine accounts for both  !
!     pressure gradients and surface tension effects in the flow. It           !
!     integrates these effects by modifying the pressure terms to include      !
!     surface tension contributions where applicable, especially in VOF        !
!     simulations.                                                             !
!   * Velocity correction: Adjusts cell-centered velocities based on various   !
!     corrections like Choi's correction and Gu's correction. This step        !
!     involves subtracting unsteady terms and body forces from the velocities. !
!   * Flux computation: Computes the mass or volume fluxes at the cell faces   !
!     by considering interpolated velocity, pressure gradients, and additional !
!     corrections. This is a crucial step in maintaining the accuracy and      !
!     stability of the flow solution.                                          !
!   * User-defined functions: Incorporates user-defined functions to allow     !
!     for specific customizations and behaviors within the interpolation       !
!     process.                                                                 !
!   * Performance monitoring: Monitors the subroutine's performance throughout !
!     its execution, contributing to the analysis and optimization of the      !
!     simulation process.                                                      !
!------------------------------------------------------------------------------!
!   Note                                                                       !
!                                                                              !
!   * Remember one important thing:                                            !
!     > sigma * kappa * dc/dx ~ dp/dx  [N / m^3]                               !
!                                                                              !
!   * Meaning also that:                                                       !
!     > sigma * kappa * (c2 - c1) ~ (p2 - p1)  [N / m^2]                       !
!                                                                              !
!   * and these terms are bundled together in this function.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Vof_Type),      target :: Vof      !! VOF object
  type(Native_Type),   target :: Nat      !! native solver object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Var_Type),    pointer :: col             ! sharp color function
  type(Matrix_Type), pointer :: A               ! pressure matrix
  type(Matrix_Type), pointer :: M               ! momentum matrix
  real                       :: a12, px_f, py_f, pz_f, fs, sigma
  real                       :: w_o, w_oo, curv_f
  integer                    :: s, c1, c2, c
  real, contiguous,  pointer :: u_c(:), v_c(:), w_c(:), v_m(:), t_m(:)
  real, contiguous,  pointer :: pst_x(:), pst_y(:), pst_z(:), pst_d(:)
  real, contiguous,  pointer :: u_f(:), v_f(:), w_f(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Rhie_And_Chow')

  call Work % Connect_Real_Cell(u_c, v_c, w_c, v_m, t_m, pst_x, pst_y, pst_z)
  call Work % Connect_Real_Face(pst_d, u_f, v_f, w_f)

  ! Take aliases
  Grid   => Flow % pnt_grid
  p      => Flow % p
  v_flux => Flow % v_flux
  col    => Vof % fun
  A      => Nat % A
  M      => Nat % M
  sigma  =  Vof % surface_tension
  call Flow % Alias_Momentum(u, v, w)

  !--------------------------------------!
  !   Store Grid % vol(c) / M % sav(c)   !
  !--------------------------------------!
  ! Units here: m^3 s / kg
  ! Remember:  M   is big in liquid, small in gas
  ! meaaning:  v_m is big in gas, small in liquid
  do c = 1, Grid % n_cells
    v_m(c) = Grid % vol(c) / M % sav(c)
  end do

  !------------------------------------------------------------------!
  !   Store pressure gradients in cells and differences over faces   !
  !------------------------------------------------------------------!

  ! Pressure gradients
  do c = 1, Grid % n_cells
    pst_x(c) = p % x(c)
    pst_y(c) = p % y(c)
    pst_z(c) = p % z(c)
  end do

  ! Pressure differences
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    pst_d(s) = p % n(c1) - p % n(c2)
  end do

  !---------------------------------------------------------------------!
  !   In case of vof simulations, treat surface tension like pressure   !
  !   The most elegat way to do it - just bundle them up in one array   !
  !---------------------------------------------------------------------!
  if(Flow % with_interface .and.  &
     Vof % surface_tension > TINY) then

    ! Surface tension gradients
    do c = 1, Grid % n_cells
      pst_x(c) = pst_x(c) - sigma * Vof % curv(c) * col % x(c)
      pst_y(c) = pst_y(c) - sigma * Vof % curv(c) * col % y(c)
      pst_z(c) = pst_z(c) - sigma * Vof % curv(c) * col % z(c)
    end do

    ! Surface tension differences
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      fs = Grid % f(s)

      ! Curvature at the face; unit: [1/m]
      curv_f = fs * Vof % curv(c1) + (1.0-fs) * Vof % curv(c2)
      pst_d(s) = pst_d(s) - sigma * curv_f * (col % n(c1) - col % n(c2))
    end do

  end if

  !--------------------------------------!
  !   Take velocities as last computed   !
  !--------------------------------------!
  do c = 1, Grid % n_cells
    u_c(c) = Flow % u % n(c)
    v_c(c) = Flow % v % n(c)
    w_c(c) = Flow % w % n(c)
  end do

  !--------------------------------------------------------------------------!
  !   Choi's correction, part 1: subtract the cell-centered unsteady terms   !
  !--------------------------------------------------------------------------!
  if(Flow % choi_correction) then

    ! To make sure that you have both o and oo values
    if(Time % Curr_Dt() - Time % First_Dt() > 3) then

      ! Weights of o and oo time step depending on the scheme used
      if(u % td_scheme == LINEAR) then
        w_o  = 1.0
        w_oo = 0.0
      else if(u % td_scheme == PARABOLIC) then
        w_o  =  2.0
        w_oo = -0.5
      end if

      do c = 1, Grid % n_cells
        ! Unit for t_m: m^3 * kg/m^3 / s * s/kg = 1
        t_m(c) = (Grid % vol(c) * Flow % density(c) / Flow % dt) / M % sav(c)

        u_c(c) = u_c(c) - (w_o * u % o(c) + w_oo * u % oo(c)) * t_m(c)
        v_c(c) = v_c(c) - (w_o * v % o(c) + w_oo * v % oo(c)) * t_m(c)
        w_c(c) = w_c(c) - (w_o * w % o(c) + w_oo * w % oo(c)) * t_m(c)
      end do

    end if

  end if

  !---------------------------------------------------------------------!
  !   Gu's correction, part 1: subtract the cell-centered body forces   !
  !---------------------------------------------------------------------!
  if(Flow % gu_correction) then

    do c = 1, Grid % n_cells

      ! Units: m^3 s/kg * kg /(m^2 s^2) = m / s
      u_c(c) = u_c(c) - v_m(c) * Flow % cell_fx(c)
      v_c(c) = v_c(c) - v_m(c) * Flow % cell_fy(c)
      w_c(c) = w_c(c) - v_m(c) * Flow % cell_fz(c)

    end do
  end if

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    fs = Grid % f(s)

    Assert(c2 > 0)

    ! Face is inside the domain
    if(c2 > 0) then

      ! Interpolate velocity 
      u_f(s) = fs * u_c(c1) + (1.0 - fs) * u_c(c2)
      v_f(s) = fs * v_c(c1) + (1.0 - fs) * v_c(c2)
      w_f(s) = fs * w_c(c1) + (1.0 - fs) * w_c(c2)

      ! Calculate coeficients for the pressure matrix
      ! Units: m * m^3 s / kg = m^4 s / kg
      a12 = A % fc(s) * (fs * v_m(c1) + (1.0-fs) * v_m(c2))
      A % val(A % pos(1,s)) = -a12
      A % val(A % pos(2,s)) = -a12
      A % val(A % dia(c1))  = A % val(A % dia(c1)) +  a12
      A % val(A % dia(c2))  = A % val(A % dia(c2)) +  a12

      ! Interpolate pressure gradients
      ! Units: kg/(m^2 s^2) * m^3 s / kg * m = m^2 / s
      ! Remember from above: v_m is big in gas, small in liquid, meaning
      ! that this interpolation is gas-biased, which is a good thing
      px_f = (       fs  * (pst_x(c1)*v_m(c1))    &
              + (1.0-fs) * (pst_x(c2)*v_m(c2)) )  &
            * Grid % dx(s)
      py_f = (       fs  * (pst_y(c1)*v_m(c1))    &
              + (1.0-fs) * (pst_y(c2)*v_m(c2)) )  &
            * Grid % dy(s)
      pz_f = (       fs  * (pst_z(c1)*v_m(c1))    &
              + (1.0-fs) * (pst_z(c2)*v_m(c2)) )  &
            * Grid % dz(s)

      ! Calculate mass or volume flux through cell face
      ! Units in lines which follow:
      ! m^3 / s = m/s * m^2
      !         + m^4 s / kg * kg / (m s^2)
      !         + m * m^2/s = m^3/s
      v_flux % n(s) = (  u_f(s) * Grid % sx(s)             &
                       + v_f(s) * Grid % sy(s)             &
                       + w_f(s) * Grid % sz(s) )           &
                       + a12 * pst_d(s)                    &
                       + A % fc(s) * (px_f + py_f + pz_f)

      !------------------------------------------------------------!
      !   Choi's correction, part 2: add flux from old time step   !
      !------------------------------------------------------------!
      if(Flow % choi_correction) then

        ! To make sure that you have both o and oo values
        if(Time % Curr_Dt() - Time % First_Dt() > 3) then

          v_flux % n(s) = v_flux % n(s)                        &
                        + (fs * t_m(c1) + (1.0-fs) * t_m(c2))  &
                          * (w_o * v_flux % o(s) + w_oo * v_flux % oo(s))
        end if
      end if  ! choi_correction

      !-------------------------------------------------------!
      !   Gu's correction, part 2: add face-centered forces   !
      !-------------------------------------------------------!
      ! Units: m^3 s/kg * kg/(m^2 s^2) * m^2 = m^3/s
      if(Flow % gu_correction) then
        v_flux % n(s) = v_flux % n(s)                          &
                      + (fs * v_m(c1) + (1.0-fs) * v_m(c2))    &
                      * (  Flow % face_fx(s) * Grid % sx(s)    &
                         + Flow % face_fy(s) * Grid % sy(s)    &
                         + Flow % face_fz(s) * Grid % sz(s) )
      end if  ! gu_correction

    end if
  end do

  call Work % Disconnect_Real_Cell(u_c, v_c, w_c, v_m, t_m, pst_x, pst_y, pst_z)
  call Work % Disconnect_Real_Face(pst_d, u_f, v_f, w_f)

  call Profiler % Stop('Rhie_And_Chow')

  end subroutine
