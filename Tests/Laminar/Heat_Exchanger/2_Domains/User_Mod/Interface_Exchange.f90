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
  real    :: t1, t2, k1, wd1, k2, wd2
  integer, parameter :: T  = 1,  &  ! store temperature as the first ...
                        K  = 2,  &  ! ... conductivity as the second ...
                        WD = 3      ! ... and wall distance as the third
!==============================================================================!

  !-------------------------------------------!
  !   Store the desired values to interface   !
  !-------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom
      call Interface_Mod_To_Buffer(inter(d1, d2),        &
                                   flow(d1) % t % n,     &
                                   flow(d2) % t % n,     &
                                   T)
      call Interface_Mod_To_Buffer(inter(d1, d2),            &
                                   flow(d1) % conductivity,  &
                                   flow(d2) % conductivity,  &
                                   K)
      call Interface_Mod_To_Buffer(inter(d1, d2),                    &
                                   flow(d1) % pnt_grid % wall_dist,  &
                                   flow(d2) % pnt_grid % wall_dist,  &
                                   WD)
    end do
  end do

  !---------------------------------------!
  !   Use the values from the interface   !
  !---------------------------------------!

  do d1 = 1, n_dom
    do d2 = 1, n_dom

      ! On the side of domain 1
      ! (ic1 and bc1 here mean cell inside and on the boundary of domain 2)
      do n1 = 1, inter(d1, d2) % n1_sub
        n   = inter(d1, d2) % face_1(n1)   ! interface index
        ic1 = inter(d1, d2) % cell_1(n1)   ! domain 1, cell inside the domain
        bc1 = inter(d1, d2) % bcel_1(n1)   ! domain 1, cell on the boundary
        t1  = flow(d1) % t % n(ic1)                 ! temperature in dom 1
        k1  = flow(d1) % conductivity(ic1)          ! conductivity in dom 1
        wd1 = flow(d1) % pnt_grid % wall_dist(ic1)  ! wall distance in dom 1
        t2  = inter(d1, d2) % phi_2(n, T)           ! temperature in dom 2
        k2  = inter(d1, d2) % phi_2(n, K)           ! conductivity in dom 2
        wd2 = inter(d1, d2) % phi_2(n, WD)          ! wall distance in dom 2
        flow(d1) % t % n(bc1) = (t1 * k1 / wd1 + t2 * k2 / wd2)  &
                              / (     k1 / wd1 +      k2 / wd2)
      end do

      ! On the side of domain 2
      ! (ic2 and bc2 here mean cell inside and on the boundary of domain 2)
      do n2 = 1, inter(d1, d2) % n2_sub
        n   = inter(d1, d2) % face_2(n2)   ! interface index
        ic2 = inter(d1, d2) % cell_2(n2)   ! domain 1, cell inside the domain
        bc2 = inter(d1, d2) % bcel_2(n2)   ! domain 1, cell on the boundary
        t2  = flow(d2) % t % n(ic2)                 ! temperature in dom 2
        k2  = flow(d2) % conductivity(ic2)          ! conductivity in dom 2
        wd2 = flow(d2) % pnt_grid % wall_dist(ic2)  ! wall distance in dom 2
        t1  = inter(d1, d2) % phi_1(n, T)           ! temperature in dom 1
        k1  = inter(d1, d2) % phi_1(n, K)           ! conductivity in dom 1
        wd1 = inter(d1, d2) % phi_1(n, WD)          ! wall distance in dom 1
        flow(d2) % t % n(bc2) = (t1 * k1 / wd1 + t2 * k2 / wd2)  &
                              / (     k1 / wd1 +      k2 / wd2)
      end do
    end do
  end do

  end subroutine
