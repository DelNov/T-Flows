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
  integer                  :: d1, d2, n1, n2, s
  integer                  :: c1i, c1b, c2i, c2b
  real                     :: t_int, cond1, cond2
!==============================================================================!

  if(n_dom < 2) return

  call Cpu_Timer_Mod_Start('Interface_Mod_Exchange')

  ! Copy values from inside cells to the interface
  do d1 = 1, n_dom
    grid1 => flow(d1) % pnt_grid
    do d2 = 1, n_dom
      grid2 => flow(d2) % pnt_grid

      ! On the side of domain 1
      do s = 1, inter(d1, d2) % n_faces
        c1i = inter(d1, d2) % glo_1(s)   ! domain 1, cell 1
        inter(d1, d2) % phi_1(s) = flow(d1) % t % n(c1i)
      end do

      ! On the side of domain 2
      do s = 1, inter(d1, d2) % n_faces
        c2i = inter(d1, d2) % glo_2(s)   ! domain 2, cell 1
        inter(d1, d2) % phi_2(s) = flow(d2) % t % n(c2i)
      end do

    end do
  end do

  ! Here you will have to exchange (global sum) of phi_1 and phi_2

  ! Copy values from interface back to boundary cells
  do d1 = 1, n_dom
    grid1 => flow(d1) % pnt_grid
    do d2 = 1, n_dom
      grid2 => flow(d2) % pnt_grid

      ! On the side of domain 1
      do s = 1, inter(d1, d2) % n_faces
        c1b = inter(d1, d2) % bnd_1(s)   ! domain 1, cell 2  (bnd cell)
        flow(d1) % t % n(c1b) = (  inter(d1, d2) % phi_1(s)  &
                                 + inter(d1, d2) % phi_2(s)) * 0.5
      end do

      ! On the side of domain 2
      do s = 1, inter(d1, d2) % n_faces
        c2b = inter(d1, d2) % bnd_2(s)   ! domain 2, cell 2  (bnd cell)
        flow(d2) % t % n(c2b) = (  inter(d1, d2) % phi_1(s)  &
                                 + inter(d1, d2) % phi_2(s)) * 0.5
      end do
    end do
  end do

  call Cpu_Timer_Mod_Stop('Interface_Mod_Exchange')

  end subroutine
