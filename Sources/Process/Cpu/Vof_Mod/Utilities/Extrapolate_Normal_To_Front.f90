!==============================================================================!
  subroutine Extrapolate_Normal_To_Front(Vof, Flow, phi, towards)
!------------------------------------------------------------------------------!
!   Extrapolates (advects) a quantity towards an interface.  Normal to the     !
!   interface acts as an advection velocity here.                              !
!                                                                              !
!   The equation we are solving reads:                                         !
!                                                                              !
!   d phi                                                                      !
!   ----- + n GRAD phi = 0            [P/m]                                    !
!   dtau                                                                       !
!                                                                              !
!   Here P is the unit of the sent parameter phi.  In reality it is [K/m],     !
!   but let's keep the comments as general as we can.  Unit for dtau is [m]    !
!   and for the GRAD [1/m].  Unit for the normal is [1].                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Field_Type),  target :: Flow
  real                      :: phi(-Flow % pnt_grid % n_bnd_cells  &
                                   :Flow % pnt_grid % n_cells)
  integer, intent(in)       :: towards
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer   :: Grid
  type(Front_Type), pointer   :: Front
  real, pointer, contiguous   :: phi_n(:), phi_o(:), adv_t(:)
  real                        :: dx, dy, dz, nx, ny, nz
  real                        :: dtau, flux_s, max_error
  integer                     :: t_iter, c, s, c1, c2
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_ITER = 120
!==============================================================================!

  call Work % Connect_Real_Cell(phi_n, phi_o, adv_t)

  ! Take aliases
  Grid  => Flow % pnt_grid
  Front => Vof % Front

  !--------------------------------------------!
  !   Find dtau in coordinate invariant mode   !
  !   Unit for dtau is [m]                     !
  !--------------------------------------------!
  dtau = HUGE
  do s = 1, Grid % n_faces
    dx = abs(Grid % dx(s))
    dy = abs(Grid % dy(s))
    dz = abs(Grid % dz(s))
    dtau = min(dtau, sqrt((dx**2 + dy**2 + dz**2)))
  end do
  call Global % Min_Real(dtau)

  !-------------------!
  !   Initial guess   !
  !-------------------!
  do c = Cells_In_Domain_And_Buffers()
    phi_n(c) = phi(c)
    phi_o(c) = phi(c)
  end do

  !--------------------------------!
  !                                !
  !   Pseudo time-step iteration   !
  !                                !
  !--------------------------------!
  do t_iter = 1, MAX_ITER

    !-------------------!
    !   Unsteady term   !
    !-------------------!
    do c = Cells_In_Domain_And_Buffers()
      phi_n(c) = phi_o(c)
    end do

    !--------------------!
    !   Advection term   !
    !--------------------!
    adv_t(:) = 0.0
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      if(Front % intersects_face(s)) then

        ! Take normal in the upwind fashion
        if(towards .eq. 1 .and. Vof % fun % n(c1) < 0.5 .or.  &
           towards .eq. 0 .and. Vof % fun % n(c1) > 0.5) then
          nx = Vof % nx(c1)
          ny = Vof % ny(c1)
          nz = Vof % nz(c1)
        else
          nx = Vof % nx(c2)
          ny = Vof % ny(c2)
          nz = Vof % nz(c2)
        end if
        Assert(abs(nx) > 0.0 .or. abs(ny) > 0.0 .or. abs(nz) > 0.0)

        ! Estimate flux. Advection velocity, in this subroutine, is normal
        ! vector to interface.  Flux at a face in T-Flows is, by definition,
        ! positive if it goes from cell c1 to c2.  Unit for flux_s is [m^2]
        flux_s = (  nx * Grid % sx(s)  &
                  + ny * Grid % sy(s)  &
                  + nz * Grid % sz(s)) / Grid % s(s)
        ! Remember that normal is positive if it points from phase 0 to 1
        ! Hence, the flux we estimated above has a valid sign if we are
        ! advecting (extrapolating) values from phase 0 to phase 1.  If
        ! That is not the case, it should change the sign.  Effectivelly,
        ! this is changing the direction of the normal ("velocity")

        if(towards .eq. 0) then
          flux_s = -flux_s
        end if

        ! The case of extrapolation from phase 0 towards phase 1.  This is
        ! upwind scheme and only the values towards which the normal/flow
        ! goes is updated, which is ensured by the max/min operators
        if(towards .eq. 1) then
          ! Value in c2 is updated when flux is positive ...
          ! ... meaning when it goes from c1 to c2 (0 --> 1)
          if(Vof % fun % n(c2) > 0.5 .and. Vof % fun % n(c1) < 0.5) then
            ! Units: [m^2] * [P] = [m^2] * [P] + [m^2] * [P]   --> good!
            adv_t(c2) = adv_t(c2)  &
                      + max(flux_s, 0.0) * (phi_n(c1) - phi_n(c2)) / Grid % d(s)
          end if

          ! Value in c1 is updated when flux is negative ...
          ! ... meaning when it goes from c2 to c1 (0 --> 1)
          if(Vof % fun % n(c1) > 0.5 .and. Vof % fun % n(c2) < 0.5) then
            ! Units: [m^2] * [P] = [m^2] * [P] + [m^2] * [P]   --> good!
            adv_t(c1) = adv_t(c1)  &
                      - min(flux_s, 0.0) * (phi_n(c2) - phi_n(c1)) / Grid % d(s)
          end if
        end if

        ! The case of extrapolation from phase 1 towards phase 0.  This is
        ! upwind scheme and only the values towards which the normal/flow
        ! goes is updated, which is ensured by the max/min operators.
        if(towards .eq. 0) then
          ! Value in c2 is updated when flux is positive ...
          ! ... meaning when goes from c1 to c2 (1 --> 0)
          if(Vof % fun % n(c2) < 0.5 .and. Vof % fun % n(c1) > 0.5) then
            ! Units: [m^2] * [P] = [m^2] * [P] + [m^2] * [P]   --> good!
            adv_t(c2) = adv_t(c2)  &
                      + max(flux_s, 0.0) * (phi_n(c1) - phi_n(c2)) / Grid % d(s)
          end if
          ! Value in c1 is updated when flux is negative ...
          ! ... meaning when it goes from c2 to c1 (1 --> 0)
          if(Vof % fun % n(c1) < 0.5 .and. Vof % fun % n(c2) > 0.5) then
            ! Units: [m^2] * [P] = [m^2] * [P] + [m^2] * [P]   --> good!
            adv_t(c1) = adv_t(c1)  &
                      - min(flux_s, 0.0) * (phi_n(c2) - phi_n(c1)) / Grid % d(s)
          end if
        end if
      end if  ! e1 > 0 .or. e2 > 0
    end do  ! through faces

    do c = Cells_In_Domain_And_Buffers()
      ! Units: [P] = [P] + [m^2] * [P] * [m] / [m^3]  --> good!
      phi_n(c) = phi_n(c) + adv_t(c) * dtau
    end do
    call Grid % Exchange_Cells_Real(phi_n)

    ! Check the error
    max_error = -HUGE
    do c = Cells_In_Domain()
      max_error = max(max_error, abs(phi_o(c) - phi_n(c)))
    end do
    call Global % Max_Real(max_error)
    if(max_error < MILI) goto 1

    ! The new becomes old
    do c = Cells_In_Domain_And_Buffers()
      phi_o(c) = phi_n(c)
    end do

  end do  ! iter

1 continue

  if(First_Proc()) then
    if(t_iter .eq. MAX_ITER) then
      PRINT *, 'Convergence not reached in the extrapolation after',  &
               t_iter, 'iterations!'
      PRINT *, 'Max error =', max_error, "Max iter", MAX_ITER
    end if
  end if

  !----------------------------------------------------------!
  !                                                          !
  !   Important: return the values to the calling function   !
  !                                                          !
  !----------------------------------------------------------!
  phi(:) = phi_n(:)

  call Work % Disconnect_Real_Cell(phi_n, phi_o, adv_t)

  end subroutine
