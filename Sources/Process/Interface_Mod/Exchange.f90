!==============================================================================!
  subroutine Interface_Mod_Exchange(inter, flow, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)     :: inter(MD, MD)
  type(Field_Type), target :: flow(MD)
  integer                  :: n_dom

  type(Grid_Type), pointer :: grid1, grid2
  integer                  :: d1, d2, n1, n2, c, s, n_tot
  real                     :: t_int, cond1, cond2
!==============================================================================!

  if(n_dom < 2) return

  call Cpu_Timer_Mod_Start('Interface_Mod_Exchange')

  !----------------------------------------------------!
  !   Copy values from inside cells to the interface   !
  !----------------------------------------------------!
  do d1 = 1, n_dom
    grid1 => flow(d1) % pnt_grid
    do d2 = 1, n_dom
      grid2 => flow(d2) % pnt_grid

      n_tot = inter(d1, d2) % n_tot

      ! Nullify values at the interface
      inter(d1, d2) % phi_1(1:n_tot) = 0.0
      inter(d1, d2) % phi_2(1:n_tot) = 0.0

      ! On the side of domain 1
      do n1 = 1, inter(d1, d2) % n1_sub
        c = inter(d1, d2) % cell_1(n1)   ! domain 1, cell inside
        s = inter(d1, d2) % face_1(n1)   ! domain 1, face
        inter(d1, d2) % phi_1(s) = flow(d1) % t % n(c)
      end do

      ! On the side of domain 2
      do n2 = 1, inter(d1, d2) % n2_sub
        c = inter(d1, d2) % cell_2(n2)   ! domain 2, cell inside
        s = inter(d1, d2) % face_2(n2)   ! domain 2, face
        inter(d1, d2) % phi_2(s) = flow(d2) % t % n(c)
      end do

      ! Here we exchange (global sum) of phi_1 and phi_2
      call Comm_Mod_Global_Sum_Real_Array(n_tot, inter(d1, d2) % phi_1(1:n_tot))
      call Comm_Mod_Global_Sum_Real_Array(n_tot, inter(d1, d2) % phi_2(1:n_tot))

    end do
  end do

  !-------------------------------------------------------!
  !   Copy values from interface back to boundary cells   !
  !-------------------------------------------------------!
  do d1 = 1, n_dom
    grid1 => flow(d1) % pnt_grid
    do d2 = 1, n_dom
      grid2 => flow(d2) % pnt_grid

      ! On the side of domain 1
      do n1 = 1, inter(d1, d2) % n1_sub
        c = inter(d1, d2) % bcel_1(n1)   ! domain 1, cell on the boundary
        s = inter(d1, d2) % face_1(n1)   ! domain 1, face
        flow(d1) % t % n(c) = (  inter(d1, d2) % phi_1(s)  &
                               + inter(d1, d2) % phi_2(s)) * 0.5
      end do

      ! On the side of domain 2
      do n2 = 1, inter(d1, d2) % n2_sub
        c = inter(d1, d2) % bcel_2(n2)   ! domain 2, cell on the boundary
        s = inter(d1, d2) % face_2(n2)   ! domain 2, face
        flow(d2) % t % n(c) = (  inter(d1, d2) % phi_1(s)  &
                               + inter(d1, d2) % phi_2(s)) * 0.5
      end do

    end do
  end do

  call Cpu_Timer_Mod_Stop('Interface_Mod_Exchange')

  end subroutine
