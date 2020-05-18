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
  real    :: u1, v1, w1
  integer, parameter :: U = 1,  &  ! store u-velocity component as first ...
                        V = 2,  &  ! ... v-velocity components as second ...
                        W = 3      ! ... and w-velocity as the third
!==============================================================================!

  !-------------------------------------------!
  !   Store the desired values to interface   !
  !-------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom
      call Interface_Mod_To_Buffer(inter(d1, d2),        &
                                   flow(d1) % u % n,     &
                                   flow(d2) % u % n,     &
                                   U)
      call Interface_Mod_To_Buffer(inter(d1, d2),        &
                                   flow(d1) % v % n,     &
                                   flow(d2) % v % n,     &
                                   V)
      call Interface_Mod_To_Buffer(inter(d1, d2),        &
                                   flow(d1) % w % n,     &
                                   flow(d2) % w % n,     &
                                   W)
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
        n   = inter(d1, d2) % face_2(n2)   ! interface index
        bc2 = inter(d1, d2) % bcel_2(n2)   ! domain 2, cell on the boundary
        u1  = inter(d1, d2) % phi_1(n, U)  ! u-velocity in domain 1
        v1  = inter(d1, d2) % phi_1(n, V)  ! u-velocity in domain 1
        w1  = inter(d1, d2) % phi_1(n, W)  ! u-velocity in domain 1
        flow(d2) % u % n(bc2) = u1
        flow(d2) % v % n(bc2) = v1
        flow(d2) % w % n(bc2) = w1
      end do
    end do
  end do

  end subroutine
