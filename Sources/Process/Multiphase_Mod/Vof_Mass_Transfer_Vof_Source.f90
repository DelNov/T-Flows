!==============================================================================!
  subroutine Multiphase_Mod_Vof_Mass_Transfer_Vof_Source(mult, b)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: b(mult % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: c, e, g, l, s, c1, c2, i_ele
  real                      :: t_x_1, t_x_2, cond_1, cond_2
!==============================================================================!

  ! Take aliases
  grid => mult % pnt_grid
  flow => mult % pnt_flow

  ! RETURN

  if(.not. flow % mass_transfer) return

  ! Distinguish between liquid and vapor
  call Multiphase_Mod_Vof_Get_Gas_And_Liquid_Phase(mult, g, l)

  call Multiphase_Mod_Vof_Calculate_Grad_Matrix_With_Front(mult)
  call Multiphase_Mod_Vof_Grad_Variable_With_Front(mult, flow % t, mult % t_sat)
  call Field_Mod_Calculate_Grad_Matrix(flow)

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(any(mult % Front % face_at_elem(1:2,s) .ne. 0)) then
      do i_ele = 1, 2
        e = mult % Front % face_at_elem(i_ele, s)
        if(e .ne. 0) then

          ! Take conductivities from each side of the interface
          cond_1 = mult % phase_cond(1)
          cond_2 = mult % phase_cond(2)
          if(mult % vof % n(c1) < 0.5) cond_1 = mult % phase_cond(2)
          if(mult % vof % n(c2) > 0.5) cond_2 = mult % phase_cond(1)

          ! Take gradients from each side of the interface
          t_x_1 = flow % t % x(c1)
          t_x_2 = flow % t % x(c2)

          IF(MATH_MOD_APPROX_REAL(GRID % YS(S), 0.0) .AND.  &
             MATH_MOD_APPROX_REAL(GRID % ZS(S), 0.0)) THEN
            WRITE(400, '(99(es12.4))')                                                           &
              t_x_1,                                                                             &
              t_x_1 * cond_1,                                                                    &
              t_x_1 * cond_1 / 2.26e+6,                                                          &
              t_x_1 * cond_1 / 2.26e+6 * (1.0/mult % phase_dens(g) - 1.0/mult % phase_dens(l)),  &
              mult % Front % elem(e) % xe
          END IF

          if(mult % Front % cell_at_elem(c1) .ne. 0) then
            mult % m_dot(c1) = -t_x_1 * cond_1 / 2.26e+6
          end if

          if(mult % Front % cell_at_elem(c2) .ne. 0) then
            mult % m_dot(c2) = -t_x_1 * cond_1 / 2.26e+6
          end if

        end if
      end do
    end if
  end do

  ! Volume source
  ! Here is the trick to get the sign correct:
  ! - if gas is 1 and liquid 2 => l-g =  1 => source > 0
  ! - if gas is 2 and liquid 1 => l-g = -1 => source < 0
  do c = 1, grid % n_cells
    e = mult % Front % cell_at_elem(c)  ! Front element
    if(e .ne. 0) then
!OK   b(c) = b(c)            &
!OK        + mult % m_dot(c) * (l-g) * mult % Front % elem(e) % area * (1.0/mult % phase_dens(g) - 1.0/mult % phase_dens(l))
      b(c) = b(c)            &
           + mult % m_dot(c) * (l-g) * mult % Front % elem(e) % area / mult % phase_dens(g)
!was OK     b(c) = b(c)            &
!was OK          - 0.0002 * mult % Front % elem(e) % area
PRINT *, 'mult % Front % elem(e) % area = ', mult % Front % elem(e) % area
    end if
  end do

  end subroutine
