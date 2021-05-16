!==============================================================================!
  subroutine Pressure_Correction(Vof, Sol)
!------------------------------------------------------------------------------!
!   Correct fluxes on pressure equation due to surface tension                 !
!                                                                              !
!   Remember one important thing:                                              !
!                                                                              !
!                   dC   dp                                                    !
!   sigma * kappa * -- ~ --  [N / m^3]                                         !
!                   dx   dx                                                    !
!                                                                              !
!   Meaning also:                                                              !
!                                                                              !
!   sigma * kappa * (C2 - C1) ~ (P2 - P1)  [N / m^2]                           !
!                                                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: pst_x => r_cell_01,  &
                      pst_y => r_cell_02,  &
                      pst_z => r_cell_03,  &
                      v_m       => r_cell_04      ! for Rhie and Chow
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Face_Type),   pointer :: v_flux
  type(Var_Type),    pointer :: col
  type(Matrix_Type), pointer :: A
  type(Matrix_Type), pointer :: M
  real, contiguous,  pointer :: b(:)
  integer                    :: c, c1, c2, s, nb, nc
  real                       :: a12, fs
  real                       :: u_fo, v_fo, w_fo
  real                       :: stens_source
  real                       :: factor2, correction, k_f, sigma
  real                       :: px_f, py_f, pz_f
!==============================================================================!

  ! Don't use this if vof is not engaged ...
  if(.not. Vof % model .eq. VOLUME_OF_FLUID) return

  ! ... or if surface tension is neglected
  if(Vof % surface_tension < TINY) return

  ! Take aliases
  grid   => Vof % pnt_grid
  flow   => Vof % pnt_flow
  col    => Vof % fun
  v_flux => flow % v_flux
  A      => Sol % A
  M      => Sol % M
  b      => Sol % b % val
  sigma  =  Vof % surface_tension
  nb     =  grid % n_bnd_cells
  nc     =  grid % n_cells
  ! Don't do what you once did in the line which follows this comment,
  ! it is plain silly as the smoothed (convoluted) variant of the vof
  ! function is used just for estimation of normals and curvatures.
  ! col    => Vof % smooth

  !--------------------------------------!
  !   Store grid % vol(c) / M % sav(c)   !
  !--------------------------------------!
  ! Units here: m^3 s / kg
  ! Remember:  M   is big in liquid, small in gas
  ! meaaning:  v_m is big in gas, small in liquid
  do c = 1, grid % n_cells
    v_m(c) = grid % vol(c) / M % sav(c)
  end do

  ! Store cell forces due to surface tension -> this is analogue to dp / dx
  ! Unit: N/m * 1/m * 1/m = N/m^3
  do c = 1, grid % n_cells
    pst_x(c) = sigma * Vof % curv(c) * col % x(c)
    pst_y(c) = sigma * Vof % curv(c) * col % y(c)
    pst_z(c) = sigma * Vof % curv(c) * col % z(c)
  end do

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    if(c2 > 0) then

      ! Interpolate pressure gradients
      ! Units: kg/(m^2 s^2) * m^3 s / kg * m = m^2 / s
      px_f = (       fs  * (pst_x(c1)*v_m(c1))    &
              + (1.0-fs) * (pst_x(c2)*v_m(c2)) )  &
            * grid % dx(s)
      py_f = (       fs  * (pst_y(c1)*v_m(c1))    &
              + (1.0-fs) * (pst_y(c2)*v_m(c2)) )  &
            * grid % dy(s)
      pz_f = (       fs  * (pst_z(c1)*v_m(c1))    &
              + (1.0-fs) * (pst_z(c2)*v_m(c2)) )  &
            * grid % dz(s)

      ! Calculate coeficients for the pressure matrix
      ! Units: m * m^3 s / kg = m^4 s / kg
      a12 = A % fc(s) * (fs * v_m(c1) + (1.0-fs) * v_m(c2))

      ! Curvature at the face; unit: [1/m]
      k_f = fs * Vof % curv(c1) + (1.0-fs) * Vof % curv(c2)

      ! Unit for stens_source: [kg/s^2 * m^4s/kg * 1/m = m^3/s]
      ! m^4 s / kg * kg / (m^2 s^2) * m = m^3/s
      stens_source = sigma * k_f * a12 * (col % n(c2)-col % n(c1))  &
                   - A % fc(s) * (px_f + py_f + pz_f)

      v_flux % n(s) = v_flux % n(s) + stens_source

      b(c1) = b(c1) - stens_source
      b(c2) = b(c2) + stens_source

    end if  ! c2 > 0

  end do

!@#  if(flow % mass_transfer) then
!@#    do c = 1, grid % n_cells
!@#      b(c) = b(c) + Vof % m_dot(c) * grid % vol(c)                    &
!@#                                    * ( 1.0 / Vof % phase_dens(1)     &
!@#                                      - 1.0 / Vof % phase_dens(2) )
!@#    end do
!@#  end if

  end subroutine
