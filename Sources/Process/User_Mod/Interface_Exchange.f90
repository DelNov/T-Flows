!==============================================================================!
  subroutine User_Mod_Interface_Exchange(inter, Flow, Turb, Vof, Swarm, n_dom)
!------------------------------------------------------------------------------!
!>  User_Mod_Interface_Exchange is a critical subroutine for simulations
!>  involving multiple domains, where it manages the exchange of data across
!>  domain interfaces.
!------------------------------------------------------------------------------!
!   Functionality (for the case of conjugate heat transfer)                    !
!                                                                              !
!   * Verification of the presence of multiple domains and heat transfer.      !
!   * Transference of variables (such as temperature, conductivity, and wall   !
!     distance) to interface buffers through Interface_Mod_To_Buffer calls.    !
!     Note that the Interface_Mod_To_Buffer also exchanges the values in the   !
!     buffers.  After this call, each domain will have values from             !
!     neighbouring domains at her disposal.                                    !
!   * Application of buffered values to impose tailored boundary conditions    !
!     within each domain, crucial for simulations with interactive domain      !
!     behavior, like conjugate heat transfer.                                  !
!------------------------------------------------------------------------------!
!   Note                                                                       !
!                                                                              !
!   * Although this example deals with the case of conjugate heat transfer     !
!     problems only, it gives a general overview on how to deal with exchange  !
!     of data at domain interfaces.  Another useful example can be found here: !
!     [root]/Tests/Laminar/Copy_Inlet/User_Mod/Interface_Exchange.f90          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)     :: inter(MD, MD)  !! parent Interface_Type
  type(Field_Type), target :: Flow(MD)       !! flows involved in simulation
  type(Turb_Type),  target :: Turb(MD)       !! turbulent fields in simulation
  type(Vof_Type),   target :: Vof(MD)        !! volume of fluid functions
  type(Swarm_Type), target :: Swarm(MD)      !! swarms of particles
  integer,      intent(in) :: n_dom          !! number of domains
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: T  = 1,  &  ! store temperature as the first ...
                        K  = 2,  &  ! ... conductivity as the second ...
                        WD = 3      ! ... and wall distance as the third.
!-----------------------------------[Locals]-----------------------------------!
  integer :: d1, d2     ! counters for domains 1 and 2
  integer :: n1, n2, n  ! local counters and number at interfaces 1 and 2
  integer :: ic1, bc1   ! internal (i) and boundary (b) cells in domain 1
  integer :: ic2, bc2   ! internal (i) and boundary (b) cells in domain 1
  real    :: t1, t2     ! temperatures in domains 1 and 2
  real    :: k1, k2     ! conductivities in domains 1 and 2
  real    :: wd1, wd2   ! wall distances in domains 1 and 2
!==============================================================================!

  ! If there is only one domain, there is nothing to exchange
  if(n_dom < 2) return

  ! If there is no heat transfer, there is nothing to exchange
  do d1 = 1, n_dom
    if(.not. Flow(d1) % heat_transfer) return
  end do

  !-------------------------------------------!
  !                                           !
  !   Store the desired values to interface   !
  !                                           !
  !-------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom

      ! Send temperature to interface
      call Interface_Mod_To_Buffer(inter(d1, d2),        &
                                   Flow(d1) % t % n,     &
                                   Flow(d2) % t % n,     &
                                   T)

      ! Send conductivities as well
      call Interface_Mod_To_Buffer(inter(d1, d2),            &
                                   Flow(d1) % conductivity,  &
                                   Flow(d2) % conductivity,  &
                                   K)

      ! Send wall distance to the interface
      call Interface_Mod_To_Buffer(inter(d1, d2),                    &
                                   Flow(d1) % pnt_grid % wall_dist,  &
                                   Flow(d2) % pnt_grid % wall_dist,  &
                                   WD)
    end do
  end do

  !------------------------------------------------------!
  !                                                      !
  !   Use the values you sent to the interface buffers   !
  !    to impose boundary conditions for each domain.    !
  !                                                      !
  !------------------------------------------------------!

  do d1 = 1, n_dom
    do d2 = 1, n_dom

      !-----------------------------!
      !   On the side of domain 1   !
      !-----------------------------!

      ! (ic1 and bc1 here mean cell inside and on the boundary of domain 2)
      do n1 = 1, inter(d1, d2) % n1_sub

        ! Fetch indexes
        n   = inter(d1, d2) % face_1(n1)   ! interface index
        ic1 = inter(d1, d2) % cell_1(n1)   ! domain 1, cell inside the domain
        bc1 = inter(d1, d2) % bcel_1(n1)   ! domain 1, cell on the boundary

        ! Fetch dependent variables from domain 1 (this domain)
        t1  = Flow(d1) % t % n(ic1)                 ! temperature in dom 1
        k1  = Flow(d1) % conductivity(ic1)          ! conductivity in dom 1
        wd1 = Flow(d1) % pnt_grid % wall_dist(ic1)  ! wall distance in dom 1

        ! Fetch values from buffers (other domain)
        t2  = inter(d1, d2) % phi_2(n, T)           ! temperature in dom 2
        k2  = inter(d1, d2) % phi_2(n, K)           ! conductivity in dom 2
        wd2 = inter(d1, d2) % phi_2(n, WD)          ! wall distance in dom 2

        !---------------------------------------------!
        !   Implementation of your model comes here   !
        !---------------------------------------------!

        ! Set temperature at the boundary of domain 1
        Flow(d1) % t % n(bc1) = (t1 * k1 / wd1 + t2 * k2 / wd2)  &
                              / (     k1 / wd1 +      k2 / wd2)

      end do

      !-----------------------------!
      !   On the side of domain 2   !
      !-----------------------------!

      ! (ic2 and bc2 here mean cell inside and on the boundary of domain 2)
      do n2 = 1, inter(d1, d2) % n2_sub

        ! Fetch indexes
        n   = inter(d1, d2) % face_2(n2)   ! interface index
        ic2 = inter(d1, d2) % cell_2(n2)   ! domain 1, cell inside the domain
        bc2 = inter(d1, d2) % bcel_2(n2)   ! domain 1, cell on the boundary

        ! Fetch dependent variables
        t2  = Flow(d2) % t % n(ic2)                 ! temperature in dom 2
        k2  = Flow(d2) % conductivity(ic2)          ! conductivity in dom 2
        wd2 = Flow(d2) % pnt_grid % wall_dist(ic2)  ! wall distance in dom 2

        ! Fetch values from buffers (other domain in this case 1)
        t1  = inter(d1, d2) % phi_1(n, T)           ! temperature in dom 1
        k1  = inter(d1, d2) % phi_1(n, K)           ! conductivity in dom 1
        wd1 = inter(d1, d2) % phi_1(n, WD)          ! wall distance in dom 1

        !---------------------------------------------!
        !   Implementation of your model comes here   !
        !---------------------------------------------!

        ! Set temperature at the boundary of domain 2
        Flow(d2) % t % n(bc2) = (t1 * k1 / wd1 + t2 * k2 / wd2)  &
                              / (     k1 / wd1 +      k2 / wd2)

      end do
    end do
  end do

  end subroutine
