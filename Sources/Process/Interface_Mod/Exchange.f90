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

  type(Grid_Type), pointer :: grid_1, grid_2
  integer                  :: d_1, d_2, n_1, n_2, s_1, s_2
  integer                  :: c_11, c_12, c_21, c_22
  real                     :: t_int, cond_1, cond_2
!==============================================================================!

  call Cpu_Timer_Mod_Start('Interface_Mod_Exchange')

  ! Write some debugging info
  do d_1 = 1, n_dom
    grid_1 => flow(d_1) % pnt_grid
    do d_2 = 1, n_dom
      grid_2 => flow(d_2) % pnt_grid
      if(inter(d_1, d_2) % n_faces > 0) then
        do n_1 = 1, inter(d_1, d_2) % n_faces
          s_1 = inter(d_1, d_2) % faces_1(n_1)     ! face in dom 1
          n_2 = inter(d_1, d_2) % close_in_2(n_1)
          s_2 = inter(d_1, d_2) % faces_2(n_2)     ! face in dom 2
          c_11 = grid_1 % faces_c(1, s_1)
          c_12 = grid_1 % faces_c(2, s_1)
          c_21 = grid_2 % faces_c(1, s_2)
          c_22 = grid_2 % faces_c(2, s_2)

          cond_1 = flow(d_1) % conductivity(c_11)
          cond_2 = flow(d_2) % conductivity(c_21)

          ! Impose interface condition (assuming delta_1 == delta_2)
          t_int = (  flow(d_1) % t % n(c_11) * cond_1    &
                   + flow(d_2) % t % n(c_21) * cond_2 )  &
                / ( cond_1 + cond_2)
          flow(d_1) % t % n(c_12) = t_int
          flow(d_2) % t % n(c_22) = t_int

!         WRITE(100,'(4F9.4)')                                 &
!           flow(d_1) % t % n(c_11), flow(d_1) % t % n(c_12),  &
!           flow(d_2) % t % n(c_21), flow(d_2) % t % n(c_22)
        end do
      end if
    end do
  end do

  call Cpu_Timer_Mod_Stop('Interface_Mod_Exchange')

  end subroutine
