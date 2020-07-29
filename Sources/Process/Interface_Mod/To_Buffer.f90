!==============================================================================!
  subroutine Interface_Mod_To_Buffer(inter, var1, var2, v)
!------------------------------------------------------------------------------!
!   Copy values of specified variables to buffer.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type) :: inter
  real                 :: var1(-inter % pnt_grid1 % n_bnd_cells :  &
                                inter % pnt_grid1 % n_cells)
  real                 :: var2(-inter % pnt_grid2 % n_bnd_cells :  &
                                inter % pnt_grid2 % n_cells)
  integer              :: v              ! variable number
!-----------------------------------[Locals]-----------------------------------!
  integer :: n1, n2, c, n, n_tot
!==============================================================================!

  !----------------------------------------------------!
  !   Copy values from inside cells to the interface   !
  !----------------------------------------------------!

  n_tot = inter % n_tot

  if(n_tot > 0) then

    ! Nullify values at the interface
    inter % phi_1(1:n_tot, v) = 0.0
    inter % phi_2(1:n_tot, v) = 0.0

    ! On the side of domain 1
    do n1 = 1, inter % n1_sub
      c = inter % cell_1(n1)   ! domain 1, cell inside
      n = inter % face_1(n1)   ! domain 1, face
      inter % phi_1(n, v) = var1(c)
    end do

    ! On the side of domain 2
    do n2 = 1, inter % n2_sub
      c = inter % cell_2(n2)   ! domain 2, cell inside
      n = inter % face_2(n2)   ! domain 2, face
      inter % phi_2(n, v) = var2(c)
    end do

    ! Here we exchange (global sum) of phi_1 and phi_2
    call Comm_Mod_Global_Sum_Real_Array(n_tot, inter % phi_1(1:n_tot,v))
    call Comm_Mod_Global_Sum_Real_Array(n_tot, inter % phi_2(1:n_tot,v))

  end if

  end subroutine
