!==============================================================================!
  subroutine User_Mod_Interface_Exchange(inter, flow, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)        :: inter(MD, MD)
  type(Field_Type),    target :: flow(MD)
  integer                     :: n_dom
!-----------------------------------[Locals]-----------------------------------!
  integer :: d1, d2, n1, n2, n, ic1, bc1, ic2, bc2
  real    :: t1, t2, k1, k2, wd1, wd2, sc1, sc2
  integer, parameter :: T  = 1,  &  ! store temperature as the first ...
                        K  = 2,  &  ! ... conductivity as the second ...
                        WD = 3,  &  ! ... and wall distance as the third ...
                        C  = 4      ! ... and scalar as the fourth
!==============================================================================!

  !-------------------------------------------!
  !                                           !
  !   Store the desired values to interface   !
  !                                           !
  !-------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom

      ! Send temperature to interface
      call Interface_Mod_To_Buffer(inter(d1, d2),        &
                                   flow(d1) % t % n,     &
                                   flow(d2) % t % n,     &
                                   T)

      ! Send conductivities as well
      call Interface_Mod_To_Buffer(inter(d1, d2),            &
                                   flow(d1) % conductivity,  &
                                   flow(d2) % conductivity,  &
                                   K)

      ! Send wall distance to the interface
      call Interface_Mod_To_Buffer(inter(d1, d2),                    &
                                   flow(d1) % pnt_grid % wall_dist,  &
                                   flow(d2) % pnt_grid % wall_dist,  &
                                   WD)

      ! Send scalar 1 as well
      call Interface_Mod_To_Buffer(inter(d1, d2),             &
                                   flow(d1) % scalar(1) % n,  &
                                   flow(d2) % scalar(1) % n,  &
                                   C)
    end do
  end do

  !---------------------------------------!
  !                                       !
  !   Use the values from the interface   !
  !                                       !
  !---------------------------------------!

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

        ! Fetch dependent variables
        t1  = flow(d1) % t % n(ic1)                 ! temperature in dom 1
        k1  = flow(d1) % conductivity(ic1)          ! conductivity in dom 1
        wd1 = flow(d1) % pnt_grid % wall_dist(ic1)  ! wall distance in dom 1
        sc1 = flow(d1) % scalar(1) % n(ic1)         ! scalar in dom 1

        t2  = inter(d1, d2) % phi_2(n, T)           ! temperature in dom 2
        k2  = inter(d1, d2) % phi_2(n, K)           ! conductivity in dom 2
        wd2 = inter(d1, d2) % phi_2(n, WD)          ! wall distance in dom 2
        sc2 = inter(d1, d2) % phi_2(n, C)           ! scalar in dom 2

        ! Set temperature at the boundary of domain 1
        flow(d1) % t % n(bc1) = (t1 * k1 / wd1 + t2 * k2 / wd2)  &
                              / (     k1 / wd1 +      k2 / wd2)

        ! Set scalar at the boundary of domain 1
        flow(d1) % scalar(1) % n(bc1) = (sc1 + sc2) * 0.5
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
        t2  = flow(d2) % t % n(ic2)                 ! temperature in dom 2
        k2  = flow(d2) % conductivity(ic2)          ! conductivity in dom 2
        wd2 = flow(d2) % pnt_grid % wall_dist(ic2)  ! wall distance in dom 2
        sc2 = flow(d2) % scalar(1) % n(ic2)         ! scalar in dom 2

        t1  = inter(d1, d2) % phi_1(n, T)           ! temperature in dom 1
        k1  = inter(d1, d2) % phi_1(n, K)           ! conductivity in dom 1
        wd1 = inter(d1, d2) % phi_1(n, WD)          ! wall distance in dom 1
        sc1 = inter(d1, d2) % phi_1(n, C)           ! scalar in dom 1

        ! Set temperature at the boundary of domain 2
        flow(d2) % t % n(bc2) = (t1 * k1 / wd1 + t2 * k2 / wd2)  &
                              / (     k1 / wd1 +      k2 / wd2)

        ! Set scalar at the boundary of domain 2
        flow(d2) % scalar(1) % n(bc2) = (sc1 + sc2) * 0.5
      end do
    end do
  end do

  end subroutine
