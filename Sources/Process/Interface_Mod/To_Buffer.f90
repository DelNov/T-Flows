!==============================================================================!
  subroutine Interface_Mod_To_Buffer(inter, grid1, grid2, var1, var2, d1, d2, v)
!------------------------------------------------------------------------------!
!   Copy values of specified variables to buffer.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type) :: inter(MD, MD)
  type(Grid_Type)      :: grid1, grid2
  real                 :: var1(-grid1 % n_bnd_cells : grid1 % n_cells)
  real                 :: var2(-grid2 % n_bnd_cells : grid2 % n_cells)
  integer              :: d1, d2
  integer              :: v  ! variable number
!-----------------------------------[Locals]-----------------------------------!
  integer :: n1, n2, c, n, n_tot
!==============================================================================!

  !----------------------------------------------------!
  !   Copy values from inside cells to the interface   !
  !----------------------------------------------------!

  n_tot = inter(d1, d2) % n_tot

  if(n_tot > 0) then

    ! Nullify values at the interface
    inter(d1, d2) % phi_1(1:n_tot, v) = 0.0
    inter(d1, d2) % phi_2(1:n_tot, v) = 0.0

    ! On the side of domain 1
    do n1 = 1, inter(d1, d2) % n1_sub
      c = inter(d1, d2) % cell_1(n1)   ! domain 1, cell inside
      n = inter(d1, d2) % face_1(n1)   ! domain 1, face
      inter(d1, d2) % phi_1(n, v) = var1(c)
    end do

    ! On the side of domain 2
    do n2 = 1, inter(d1, d2) % n2_sub
      c = inter(d1, d2) % cell_2(n2)   ! domain 2, cell inside
      n = inter(d1, d2) % face_2(n2)   ! domain 2, face
      inter(d1, d2) % phi_2(n, v) = var2(c)
    end do

    ! Here we exchange (global sum) of phi_1 and phi_2
    call Comm_Mod_Global_Sum_Real_Array(n_tot, inter(d1,d2) % phi_1(1:n_tot,v))
    call Comm_Mod_Global_Sum_Real_Array(n_tot, inter(d1,d2) % phi_2(1:n_tot,v))

  end if

  end subroutine
