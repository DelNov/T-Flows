!==============================================================================!
  subroutine User_Mod_Interface_Exchange(inter, Flow, Turb, Vof, Swarm, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)        :: inter(MD, MD)
  type(Field_Type),    target :: Flow(MD)
  type(Turb_Type),     target :: Turb(MD)
  type(Vof_Type),      target :: Vof(MD)
  type(Swarm_Type),    target :: Swarm(MD)
  integer, intent(in)         :: n_dom
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: U   = 1,  &  ! store u-velocity component as first ...
                        V   = 2,  &  ! ... v-velocity components as second ...
                        W   = 3,  &  ! ... and w-velocity as the third ...
                        KIN = 4,  &  ! ... kinetic energy as the fourth ...
                        EPS = 5      ! ... and its dissipation as the fifth
!-----------------------------------[Locals]-----------------------------------!
  integer :: d1, d2, n1, n2, n, ic1, bc1, ic2, bc2
  real    :: u1, v1, w1, k1, e1
!==============================================================================!

  !-------------------------------------------!
  !   Store the desired values to interface   !
  !-------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom
      call Interface_Mod_Exchange(inter(d1, d2),        &
                                  Flow(d1) % u % n,     &
                                  Flow(d2) % u % n,     &
                                  U)
      call Interface_Mod_Exchange(inter(d1, d2),        &
                                  Flow(d1) % v % n,     &
                                  Flow(d2) % v % n,     &
                                  V)
      call Interface_Mod_Exchange(inter(d1, d2),        &
                                  Flow(d1) % w % n,     &
                                  Flow(d2) % w % n,     &
                                  W)
      call Interface_Mod_Exchange(inter(d1, d2),        &
                                  Turb(d1) % kin % n,   &
                                  Turb(d2) % kin % n,   &
                                  KIN)
      call Interface_Mod_Exchange(inter(d1, d2),        &
                                  Turb(d1) % eps % n,   &
                                  Turb(d2) % eps % n,   &
                                  EPS)
    end do
  end do

  !---------------------------------------!
  !   Use the values from the interface   !
  !---------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom

      ! Consider only domain 2 for the copy thing
      ! (Because domain 2 is the one whose boundary values will change)
      do n2 = 1, inter(d1, d2) % n2_sub
        n   = inter(d1, d2) % face_2(n2)     ! interface index
        bc2 = inter(d1, d2) % bcel_2(n2)     ! domain 2, cell on the boundary
        u1  = inter(d1, d2) % phi_1(n, U)    ! u-velocity in domain 1
        v1  = inter(d1, d2) % phi_1(n, V)    ! v-velocity in domain 1
        w1  = inter(d1, d2) % phi_1(n, W)    ! w-velocity in domain 1
        k1  = inter(d1, d2) % phi_1(n, KIN)  ! kinetic energy in domain 1
        e1  = inter(d1, d2) % phi_1(n, EPS)  ! dissipation in domain 1
        Flow(d2) % u   % n(bc2) = u1
        Flow(d2) % v   % n(bc2) = v1
        Flow(d2) % w   % n(bc2) = w1
        Turb(d2) % kin % n(bc2) = k1
        Turb(d2) % eps % n(bc2) = e1
      end do
    end do
  end do

  end subroutine
