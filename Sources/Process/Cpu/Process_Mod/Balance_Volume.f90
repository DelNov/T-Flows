!==============================================================================!
  subroutine Balance_Volume(Process, Flow, Vof)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for modifying volume fluxes and velocities
!>  at outflows to ensure volume conservation. It is crucial for maintaining
!>  the accuracy of simulations, especially in the context of pressure
!>  calculations where strict volume conservation is necessary to achieve
!>  convergence of the pressure linear solver.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Sets up necessary pointers and variables, including      !
!     grid, flow field, and volume fluxes. Establishes aliases for momentum    !
!     variables.                                                               !
!   * Boundary flux update: Updates fluxes at boundaries based on the latest   !
!     velocities. These updates may not inherently conserve volume,            !
!     necessitating further adjustments.                                       !
!   * Phase change volume adjustment: If mass transfer is involved (like in    !
!     VOF simulations), adjusts the volume source term to account for phase    !
!     change induced volume changes.                                           !
!   * Outflow volume calculation: Calculates the total volume flux and area    !
!     of outflow boundaries, which are essential for balancing the total       !
!     volume in and out of the computational domain.                           !
!   * Volume flux balance: Computes all boundary volume fluxes (both in and    !
!     out) and adjusts velocities and fluxes at outflow boundaries to balance  !
!     the total volume flux, ensuring volume conservation.                     !
!   * Velocity and flux correction: Applies corrections to velocities and      !
!     fluxes at outflow boundaries based on calculated factors to maintain     !
!     volume balance across the computational domain.                          !
!   * Buffer refresh: Updates the buffers for velocity components due to       !
!     modifications in velocity fields during the volume balancing process.    !
!   * Performance monitoring: The subroutine monitors its execution            !
!     performance to optimize simulation processes and identify potential      !
!     computational bottlenecks.                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process  !! parent class
  type(Field_Type),    target :: Flow     !! flow object
  type(Vof_Type),      target :: Vof      !! VOF object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w
  type(Face_Type), pointer :: v_flux
  integer                  :: s, c2, reg
  real                     :: fac, vol_outflow, area_outflow
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Balance_Volume')

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  v_flux => Flow % v_flux
  call Flow % Alias_Momentum(u, v, w)

  !--------------------------------------------------------!
  !                                                        !
  !   Update fluxes at boundaries with latest velocities   !
  !     (These might not obey the volume conservation)     !
  !                                                        !
  !--------------------------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. SYMMETRY) then
      do s = Faces_In_Region(reg)
        v_flux % n(s) = 0.0
      end do
    else
      do s = Faces_In_Region(reg)
        c2 = Grid % faces_c(2,s)

        v_flux % n(s) = ( u % n(c2) * Grid % sx(s)     &
                        + v % n(c2) * Grid % sy(s)     &
                        + w % n(c2) * Grid % sz(s) )
      end do  ! faces
    end if    ! boundary condition type
  end do      ! regions

  !--------------------------------------!
  !                                      !
  !   Added volume due to phase change   !
  !                                      !
  !--------------------------------------!
  if(Flow % mass_transfer_model .ne. NO_MASS_TRANSFER) then
    call Vof % Mass_Transfer_Added_Volume(bulk % vol_src)
  end if

  !----------------------------------------------------------!
  !                                                          !
  !   Check if anything comes from real outflow boundaries   !
  !   at all.  For example, at the beginning of simulation   !
  !                                                          !
  !----------------------------------------------------------!
  vol_outflow  = 0.0
  area_outflow = 0.0
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. PRESSURE .or.  &
       Grid % region % type(reg) .eq. OUTFLOW  .or.  &
       Grid % region % type(reg) .eq. CONVECT) then
      do s = Faces_In_Region(reg)

        ! Integrate volume flux only where something does come out
        if(v_flux % n(s) > 0.0) then
          vol_outflow  = vol_outflow  + v_flux % n(s)
        end if

        ! Integrate area everywhere where outflow is specified
        area_outflow = area_outflow + Grid % s(s)

      end do
    end if
  end do
  call Global % Sum_Real(vol_outflow)
  call Global % Sum_Real(area_outflow)

  !------------------------------------------------------!
  !                                                      !
  !   Calculate all boundary volume fluxes; in and out   !
  !                                                      !
  !------------------------------------------------------!
  bulk % vol_in   = 0.0
  bulk % vol_out  = 0.0
  bulk % area_in  = 0.0
  bulk % area_out = 0.0

  !---------------------------------------------------------------!
  !   This works better if domain features pressure outflow ...   !
  !---------------------------------------------------------------!
  if(Flow % has_pressure) then

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. INFLOW    .or.  &
         Grid % region % type(reg) .eq. OUTFLOW   .or.  &
         Grid % region % type(reg) .eq. CONVECT   .or.  &
         Grid % region % type(reg) .eq. PRESSURE) then
        do s = Faces_In_Region(reg)
          if(v_flux % n(s) > 0.0) then
            bulk % vol_out  = bulk % vol_out  + v_flux % n(s)
            bulk % area_out = bulk % area_out + Grid % s(s)
          else
            bulk % vol_in  = bulk % vol_in  - v_flux % n(s)
            bulk % area_in = bulk % area_in + Grid % s(s)
          end if
        end do    ! faces
      end if      ! boundary conditions
    end do        ! regions

  !-----------------------------------------------------------------------!
  !   ... but this variant works better for all other types of outflows   !
  !-----------------------------------------------------------------------!
  else
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. INFLOW) then
        do s = Faces_In_Region(reg)
          bulk % vol_in  = bulk % vol_in  - v_flux % n(s)
          bulk % area_in = bulk % area_in + Grid % s(s)
        end do
      end if    ! boundary condition
      if(Grid % region % type(reg) .eq. OUTFLOW   .or.  &
         Grid % region % type(reg) .eq. CONVECT   .or.  &
         Grid % region % type(reg) .eq. PRESSURE) then
        do s = Faces_In_Region(reg)
          bulk % vol_out  = bulk % vol_out  + v_flux % n(s)
          bulk % area_out = bulk % area_out + Grid % s(s)
        end do
      end if    ! boundary condition
    end do      ! regions

  end if  ! flow has pressure

  call Global % Sum_Reals(bulk % vol_in,   &
                          bulk % vol_out,  &
                          bulk % area_in,  &
                          bulk % area_out)

  ! Avoid divisions by zero for the cases without any fluid motion
  fac = 1.0
  if(bulk % vol_out .gt. FEMTO) fac = bulk % vol_in / (bulk % vol_out)
  bulk % area_out = max(bulk % area_out, FEMTO)

  !----------------------------------------------------------------------!
  !                                                                      !
  !   Correct velocities and volume fluxes at outlet to balance volume   !
  !    (Here you should distinguish two cases: if anything was coming    !
  !     out of real outflows, like in the first time step, or not.)      !
  !                                                                      !
  !----------------------------------------------------------------------!

  !---------------------------------------------------!
  !   Nothing comes from real outflows, no point in   !
  !   scaling with "fac", just add bulk corrections   !
  !---------------------------------------------------!
  if(Math % Approx_Real(vol_outflow, 0.0, tol=FEMTO)) then

    bulk % vol_out = 0.0

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW   .or.  &
         Grid % region % type(reg) .eq. CONVECT   .or.  &
         Grid % region % type(reg) .eq. PRESSURE) then
        do s = Faces_In_Region(reg)
          c2 = Grid % faces_c(2,s)

          ! Update velocity components ...
          u % n(c2) = (bulk % vol_in + bulk % vol_src) / area_outflow  &
                    * Grid % sx(s) / Grid % s(s)
          v % n(c2) = (bulk % vol_in + bulk % vol_src) / area_outflow  &
                    * Grid % sy(s) / Grid % s(s)
          w % n(c2) = (bulk % vol_in + bulk % vol_src) / area_outflow  &
                    * Grid % sz(s) / Grid % s(s)

          ! ... volume flux itself ...
          v_flux % n(s) = u % n(c2)*Grid % sx(s)    &
                        + v % n(c2)*Grid % sy(s)    &
                        + w % n(c2)*Grid % sz(s)

          ! ... and bulk volume out
          bulk % vol_out = bulk % vol_out + v_flux % n(s)
        end do  ! faces
      end if    ! boundary condition
    end do      ! region

    ! Holy mackrele: summ it up over all processors
    call Global % Sum_Real(bulk % vol_out)

  !-------------------------------------------------!
  !   Something is coming out from real outflows,   !
  !   presumably not in the first iteration of      !
  !   the first time step, do scaling with "fac"    !
  !-------------------------------------------------!
  else

    ! You should compute the "fac" now ...
    fac = bulk % vol_in / (bulk % vol_out + TINY)

    ! ... and correct all velocities
    bulk % vol_out = 0.0

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW   .or.  &
         Grid % region % type(reg) .eq. CONVECT   .or.  &
         Grid % region % type(reg) .eq. PRESSURE) then
        do s = Faces_In_Region(reg)
          c2 = Grid % faces_c(2,s)

          ! Update velocity components ...
          u % n(c2) = u % n(c2) * fac + bulk % vol_src / bulk % area_out  &
                    * Grid % sx(s) / Grid % s(s)
          v % n(c2) = v % n(c2) * fac + bulk % vol_src / bulk % area_out  &
                    * Grid % sy(s) / Grid % s(s)
          w % n(c2) = w % n(c2) * fac + bulk % vol_src / bulk % area_out  &
                    * Grid % sz(s) / Grid % s(s)

          ! ... volume flux itself ...
          v_flux % n(s) = u % n(c2)*Grid % sx(s)    &
                        + v % n(c2)*Grid % sy(s)    &
                        + w % n(c2)*Grid % sz(s)

          ! ... and bulk volume out
          bulk % vol_out = bulk % vol_out + v_flux % n(s)
        end do  ! faces
      end if    ! boundary condition
    end do      ! regions

    ! Holy mackrele: summ it up over all processors
    call Global % Sum_Real(bulk % vol_out)  ! not checked

  end if

  ! You were modifying the velocity components -> refresh their buffers
  call Grid % Exchange_Cells_Real(u % n)
  call Grid % Exchange_Cells_Real(v % n)
  call Grid % Exchange_Cells_Real(w % n)

  call Profiler % Stop('Balance_Volume')

  end subroutine
