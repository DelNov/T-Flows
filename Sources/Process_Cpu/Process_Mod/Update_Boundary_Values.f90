!==============================================================================!
  subroutine Update_Boundary_Values(Process, Flow, Turb, Vof, update)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to update variables on the boundaries (boundary
!>  cells) where needed.  It does not update volume fluxes at boundaries -
!>  which is good and you should keep it like that. Keep in minda that fluxes
!>  inside are calculated inside Compute_Pressure, and fluxes at boundaries
!>  are updated in Balance_Volume.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Sets up pointers and variables, including grid, flow     !
!     field, and turbulence models. Establishes aliases for various variables  !
!     such as momentum, energy, and turbulence quantities.                     !
!   * Update control: Depending on the 'update' parameter, selectively updates !
!     boundary values for momentum, turbulence, multiphase variables, energy,  !
!     or scalars.                                                              !
!   * Momentum update: For momentum updates, extrapolates boundary values for  !
!     velocity components (u, v, w) at outflow, pressure, or symmetry regions. !
!   * Multiphase update: If involved, updates the boundary values for the VOF  !
!     function at relevant boundaries.                                         !
!   * Turbulence update: Updates turbulence quantities like kinetic energy,    !
!     dissipation, and stress tensors at boundary faces, depending on the      !
!     turbulence model used (k-epsilon, Spalart-Allmaras, RSM, etc.).          !
!   * Energy update: For heat transfer simulations, updates temperature or     !
!     heat flux at boundaries based on wall conditions and specified           !
!     turbulence models. It includes both the 'old' and 'new' ways of handling !
!     wall temperature/flux calculations.                                      !
!   * Scalars update: For simulations with scalar transport, updates scalar    !
!     variables at boundaries considering wall conditions and diffusivity.     !
!   * Error checking: Ensures the 'update' parameter is valid and throws an    !
!     error if an invalid parameter is passed.                                 !
!   * Performance monitoring: Monitors the execution of the subroutine for     !
!     optimization and efficiency.                                             !
!------------------------------------------------------------------------------!
!   Note                                                                       !
!                                                                              !
!   * Ideally, this function should be called before each calculation of       !
!     gradients, but it is hardly so.  Maybe one could even think of calling   !
!     it from calculation of variables and updating those which match name     !
!     of the variable.                                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Turb_Type),     target :: Turb     !! turbulence object
  type(Vof_Type),      target :: Vof      !! VOF object
  character(*)                :: update   !! character switch to control
                                          !! which variables to update
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi, fun
  type(Var_Type),  pointer :: kin, eps, zeta, f22, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  integer                  :: c0, c1, c2, i_fac, s, s1, sc, reg
  real                     :: kin_vis, u_tau, dt_dn
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Update_Boundary_Values')

  ! Take aliases
  Grid => Flow % pnt_grid
  vis  => Turb % vis
  fun  => Vof % fun
  call Flow % Alias_Momentum    (u, v, w)
  call Flow % Alias_Energy      (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_T2          (t2)

  call String % To_Upper_Case(update)

  if(update .ne. 'MOMENTUM'   .and.  &
     update .ne. 'TURBULENCE' .and.  &
     update .ne. 'VOF'        .and.  &
     update .ne. 'ENERGY'     .and.  &
     update .ne. 'SCALARS'    .and.  &
     update .ne. 'ALL') then
    call Message % Error(72,                                              &
                         'Invalid parameter in function call. Exiting!',  &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  !--------------!
  !              !
  !   Momentum   !
  !              !
  !--------------!
  if( (update .eq. 'MOMENTUM' .or. update .eq. 'ALL') ) then

    ! On the boundary perform the extrapolation
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW  .or.  &
         Grid % region % type(reg) .eq. PRESSURE .or.  &
         Grid % region % type(reg) .eq. SYMMETRY) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          u % n(c2) = u % n(c1)
          v % n(c2) = v % n(c1)
          w % n(c2) = w % n(c1)
        end do  ! faces
      end if    ! boundary condition
    end do      ! region

  end if  ! update momentum

  !----------------!
  !                !
  !   Multiphase   !
  !                !
  !----------------!
  if( (update .eq. 'VOF' .or. update .eq. 'ALL') .and.  &
      Flow % with_interface ) then

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2 < 0) then

        if( Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW  .or.    &
            Grid % Bnd_Cond_Type(c2) .eq. PRESSURE .or.    &
            Grid % Bnd_Cond_Type(c2) .eq. CONVECT  .or.    &
            Grid % Bnd_Cond_Type(c2) .eq. SYMMETRY ) then
          fun % n(c2) = fun % n(c1)
        end if
      end if ! c2 < 0
    end do

  end if  ! update multiphase

  !----------------!
  !                !
  !   Turbulence   !
  !                !
  !----------------!
  if(update .eq. 'TURBULENCE' .or. update .eq. 'ALL') then

    !----------------------!
    !   K-epsilon-zeta-f   !
    !----------------------!
    if(Turb % model .eq. K_EPS_ZETA_F .or.  &
       Turb % model .eq. HYBRID_LES_RANS) then

      do reg = Boundary_Regions()

        ! Regions outflow, pressure or symmetry
        if(Grid % region % type(reg) .eq. OUTFLOW  .or.  &
           Grid % region % type(reg) .eq. PRESSURE .or.  &
           Grid % region % type(reg) .eq. SYMMETRY) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            kin  % n(c2) = kin  % n(c1)
            eps  % n(c2) = eps  % n(c1)
            zeta % n(c2) = zeta % n(c1)
            f22  % n(c2) = f22  % n(c1)
            if(Flow % heat_transfer) then
              t2  % n(c2) = t2  % n(c1)
            end if
          end do  ! faces
        end if    ! boundary condition
      end do      ! regions
    end if        ! turbulence model

    !---------------!
    !   K-epsilon   !
    !---------------!
    if(Turb % model .eq. K_EPS) then

      do reg = Boundary_Regions()

        ! Regions outflow, pressure or symmetry
        if(Grid % region % type(reg) .eq. OUTFLOW  .or.  &
           Grid % region % type(reg) .eq. PRESSURE .or.  &
           Grid % region % type(reg) .eq. SYMMETRY) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            kin % n(c2) = kin % n(c1)
            eps % n(c2) = eps % n(c1)
            if(Flow % heat_transfer) then
              t2  % n(c2) = t2  % n(c1)
            end if
          end do  ! faces
        end if    ! boundary condition
      end do      ! regions
    end if        ! turbulence model

    !----------------------!
    !   Spalart-Allmaras   !
    !----------------------!
    if(Turb % model .eq. SPALART_ALLMARAS .or.  &
       Turb % model .eq. DES_SPALART) then

      do reg = Boundary_Regions()

        ! Regions outflow, pressure or symmetry
        if(Grid % region % type(reg) .eq. OUTFLOW  .or.  &
           Grid % region % type(reg) .eq. PRESSURE .or.  &
           Grid % region % type(reg) .eq. SYMMETRY) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            vis % n(c2) = vis % n(c1)
          end do  ! faces
        end if    ! boundary condition
      end do      ! regions
    end if        ! turbulence model

    !----------------------------!
    !   Reynolds stress models   !
    !----------------------------!
    if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

      do reg = Boundary_Regions()

        ! Regions at solid walls
        if(Grid % region % type(reg) .eq. WALL .or.  &
           Grid % region % type(reg) .eq. WALLFL) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            uu  % n(c2) = 0.0
            vv  % n(c2) = 0.0
            ww  % n(c2) = 0.0
            uv  % n(c2) = 0.0
            uw  % n(c2) = 0.0
            vw  % n(c2) = 0.0
            kin % n(c2) = 0.0
            kin_vis = Flow % viscosity(c1) / Flow % density(c1)
            u_tau = kin_vis * sqrt(u % n(c1)**2 + v % n(c1)**2 + w % n(c1)**2)  &
                  / Grid % wall_dist(c1)
            Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(u_tau,            &
                                                     Grid % wall_dist(c1),  &
                                                     kin_vis,               &
                                                     0.0)
            if(Turb % model .eq. RSM_MANCEAU_HANJALIC) f22 % n(c2) = 0.0
          end do  ! faces

        ! Regions outflow, pressure or symmetry
        else if(Grid % region % type(reg) .eq. OUTFLOW  .or.  &
                Grid % region % type(reg) .eq. PRESSURE .or.  &
                Grid % region % type(reg) .eq. SYMMETRY) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            uu  % n(c2) = uu  % n(c1)
            vv  % n(c2) = vv  % n(c1)
            ww  % n(c2) = ww  % n(c1)
            uv  % n(c2) = uv  % n(c1)
            uw  % n(c2) = uw  % n(c1)
            vw  % n(c2) = vw  % n(c1)
            kin % n(c2) = kin % n(c1)
            eps % n(c2) = eps % n(c1)
            if(Turb % model .eq. RSM_MANCEAU_HANJALIC)  &
              f22 % n(c2) = f22 % n(c1)
          end do  ! faces
        end if    ! boundary condition
      end do      ! regions
    end if        ! turbulence model

  end if  ! update turbulence

  !-------------------!
  !                   !
  !   Heat transfer   !  =--> the old way
  !                   !
  !-------------------!
  if( (update .eq. 'ENERGY' .or. update .eq. 'ALL') .and.  &
      Flow % heat_transfer                          .and.  &
      .not. Flow % exp_temp_wall) then

    ! Initialize variables for global heat transfer
    Flow % heat        = 0.0     ! [W]
    Flow % heat_flux   = 0.0     ! [W/m^2]
    Flow % heated_area = 0.0     ! [m^2]

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(Cell_In_This_Proc(c1) .and. c2 < 0) then

        ! Wall temperature or heat fluxes for k-eps-zeta-f
        ! and high-re k-eps models. 
        if(Turb % model .eq. K_EPS_ZETA_F    .or.  &
           Turb % model .eq. HYBRID_LES_RANS .or.  &
           Turb % model .eq. LES_DYNAMIC     .or.  &
           Turb % model .eq. LES_WALE        .or.  &
           Turb % model .eq. K_EPS) then
          if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * Grid % wall_dist(c1)  &
                      / (Turb % con_w(c1) + TINY)
          else if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * Turb % con_w(c1)  &
                      / Grid % wall_dist(c1)
          end if

        ! Wall temperature or heat fluxes for other trubulence models
        else
          if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALLFL) then
            t % n(c2) = t % n(c1) + t % q(c2) * Grid % wall_dist(c1)  &
                      / Flow % conductivity(c1)
          else if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALL) then
            t % q(c2) = ( t % n(c2) - t % n(c1) ) * Flow % conductivity(c1)  &
                      / Grid % wall_dist(c1)
          end if

        end if

        ! Integrate heat and heated area
        Flow % heat = Flow % heat + t % q(c2) * Grid % s(s)
        if(abs(t % q(c2)) > TINY) then
          Flow % heated_area = Flow % heated_area + Grid % s(s)
        end if

        if( Var_Mod_Bnd_Cond_Type(t,c2) .eq. OUTFLOW .or.     &
            Var_Mod_Bnd_Cond_Type(t,c2) .eq. PRESSURE .or.    &
            Var_Mod_Bnd_Cond_Type(t,c2) .eq. SYMMETRY ) then
          t % n(c2) = t % n(c1)
        end if

      end if ! c2 < 0
    end do ! s = 1, Grid % n_faces

    !-----------------------------------------------!
    !   Integrate (summ) heated area, and heat up   !
    !-----------------------------------------------!
    call Global % Sum_Real(Flow % heat)
    call Global % Sum_Real(Flow % heated_area)
    Flow % heat_flux = Flow % heat / max(Flow % heated_area, TINY)

  end if  ! update energy and heat transfer

  !-------------------!
  !                   !
  !   Heat transfer   !  =--> the new way
  !                   !
  !-------------------!
  if( (update .eq. 'ENERGY' .or. update .eq. 'ALL') .and.  &
      Flow % heat_transfer                          .and.  &
      Flow % exp_temp_wall) then

    ! Initialize variables for global heat transfer
    Flow % heat        = 0.0     ! [W]
    Flow % heat_flux   = 0.0     ! [W/m^2]
    Flow % heated_area = 0.0     ! [m^2]

    ! Find "off the wall" cell c0
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      if(Cell_In_This_Proc(c1) .and. c2 < 0) then
        if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALL .or.  &
           Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALLFL) then
          do i_fac = 1, Grid % cells_n_faces(c1)
            s1 = Grid % cells_f(i_fac, c1)  ! side around c1
            if(s1 .ne. s) then

              ! Find the cell on the side opposite of wall cell c2
              c0 = Grid % faces_c(1,s1) + Grid % faces_c(2,s1) - c1

              ! Use wall distace criterion to tell if this is proper cell
              if(Grid % wall_dist(c0) > 1.25 * Grid % wall_dist(c1)) goto 1
            end if
          end do
        end if
      end if

      ! At this point, c0 is known
1     continue

      ! On the boundary perform the extrapolation
      if(Cell_In_This_Proc(c1) .and. c2 < 0) then

        ! In the "new way" the extrapolation is
        ! independent from turbulence model
        if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALLFL) then

          dt_dn = -t % q(c2) / Flow % conductivity(c1)

          ! Compute t % n(c2) by exponential fit
          call Math % Fit_Exp_Derivative_And_Two_Points(  &
                      dt_dn,                              &
                      Grid % wall_dist(c2), t % n(c2),    &
                      Grid % wall_dist(c1), t % n(c1),    &
                      Grid % wall_dist(c0), t % n(c0))

        else if(Var_Mod_Bnd_Cond_Type(t,c2) .eq. WALL) then

          ! Compute dt/dn at the wall by exponential fit
          call Math % Fit_Exp_Three_Points(               &
                      dt_dn,                              &
                      Grid % wall_dist(c2), t % n(c2),    &
                      Grid % wall_dist(c1), t % n(c1),    &
                      Grid % wall_dist(c0), t % n(c0))

          t % q(c2) = -dt_dn * Flow % conductivity(c1)

        end if

        ! Integrate heat and heated area
        Flow % heat = Flow % heat + t % q(c2) * Grid % s(s)
        if(abs(t % q(c2)) > TINY) then
          Flow % heated_area = Flow % heated_area + Grid % s(s)
        end if

        if( Var_Mod_Bnd_Cond_Type(t,c2) .eq. OUTFLOW .or.     &
            Var_Mod_Bnd_Cond_Type(t,c2) .eq. PRESSURE .or.    &
            Var_Mod_Bnd_Cond_Type(t,c2) .eq. SYMMETRY ) then
          t % n(c2) = t % n(c1)
        end if

      end if ! c2 < 0
    end do ! s = 1, Grid % n_faces

    !-----------------------------------------------!
    !   Integrate (summ) heated area, and heat up   !
    !-----------------------------------------------!
    call Global % Sum_Real(Flow % heat)
    call Global % Sum_Real(Flow % heated_area)
    Flow % heat_flux = Flow % heat / max(Flow % heated_area, TINY)

  end if  ! update energy and heat transfer

  !-------------!
  !             !
  !   Scalars   !
  !             !
  !-------------!
  if( (update .eq. 'SCALARS' .or. update .eq. 'ALL') ) then

    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)

      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        ! On the boundary perform the extrapolation
        if (c2 < 0) then

          ! Wall temperature or heat fluxes for k-eps-zeta-f
          ! and high-re k-eps models. 
          if(Turb % model .eq. K_EPS_ZETA_F    .or.  &
             Turb % model .eq. HYBRID_LES_RANS .or.  &
             Turb % model .eq. K_EPS) then

            if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
              phi % n(c2) = phi % n(c1) + phi % q(c2) * Grid % wall_dist(c1)  &
                        / (Turb % diff_w(c1) + TINY)
            else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL) then
              phi % q(c2) = ( phi % n(c2) - phi % n(c1) ) * Turb % diff_w(c1)  &
                        / Grid % wall_dist(c1)
            end if

          ! Scalar boundary for other trubulence models
          else
            if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
              phi % n(c2) = phi % n(c1) + phi % q(c2) * Grid % wall_dist(c1)  &
                          / Flow % diffusivity
            else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL) then
              phi % q(c2) = (phi % n(c2) - phi % n(c1)) * Flow % diffusivity &
                          / Grid % wall_dist(c1)
            end if ! WALL or WALLFL
          end if ! Turb. models

          if( Var_Mod_Bnd_Cond_Type(phi,c2) .eq. OUTFLOW .or.     &
              Var_Mod_Bnd_Cond_Type(phi,c2) .eq. PRESSURE .or.    &
              Var_Mod_Bnd_Cond_Type(phi,c2) .eq. SYMMETRY ) then
            phi % n(c2) = phi % n(c1)
          end if

        end if ! c2 < 0
      end do ! s = 1, Grid % n_faces
    end do ! sc = 1, Flow % n_scalars

  end if  ! update_scalars

  call Profiler % Stop('Update_Boundary_Values')

  end subroutine
