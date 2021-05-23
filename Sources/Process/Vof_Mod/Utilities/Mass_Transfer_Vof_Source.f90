!==============================================================================!
  subroutine Mass_Transfer_Vof_Source(Vof, b)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real                    :: b(Vof % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: Flow
  integer                   :: c, e, g, l, s, c1, c2, i_ele
  real                      :: t_x_1, t_x_2, cond_1, cond_2
!==============================================================================!

  ! Take aliases
  grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  ! If not a problem with mass transfer, get out of here
  if(.not. Flow % mass_transfer) return

  ! Initialize mass transfer term
  Vof % m_dot(:) = 0.0

  ! Distinguish between liquid and vapor
  call Vof % Get_Gas_And_Liquid_Phase(g, l)

  !------------------------------------------------!
  !   Compute gradients of temperature, imposing   !
  !    saturation temperature at the interface     !
  !------------------------------------------------!
  call Vof % Calculate_Grad_Matrix_With_Front()
  call Vof % Grad_Variable_With_Front(Flow % t, Vof % t_sat)
  call Flow % Calculate_Grad_Matrix()

  !----------------------------------------!
  !   Compute heat flux at the interface   !
  !----------------------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(any(Vof % Front % face_at_elem(1:2,s) .ne. 0)) then
      do i_ele = 1, 2
        e = Vof % Front % face_at_elem(i_ele, s)
        if(e .ne. 0) then

          ! Take conductivities from each side of the interface
          cond_1 = Vof % phase_cond(1)
          cond_2 = Vof % phase_cond(2)
          if(Vof % fun % n(c1) < 0.5) cond_1 = Vof % phase_cond(2)
          if(Vof % fun % n(c2) > 0.5) cond_2 = Vof % phase_cond(1)

          ! Take gradients from each side of the interface
          t_x_1 = Flow % t % x(c1)
          t_x_2 = Flow % t % x(c2)

          ! WRITE DOWN STEFAN'S SOLUTION
          IF(MATH % APPROX_REAL(GRID % YS(S), 0.0) .AND.  &
             MATH % APPROX_REAL(GRID % ZS(S), 0.0)) THEN
            WRITE(400, '(99(es12.4))')                                  &
              T_X_1,                                                    &
              T_X_1 * COND_1,                                           &
              T_X_1 * COND_1 / 2.26E+6,                                 &
              T_X_1 * COND_1 / 2.26E+6 * (  1.0/VOF % PHASE_DENS(G)     &
                                          - 1.0/VOF % PHASE_DENS(L) ),  &
              VOF % FRONT % ELEM(E) % XE
          END IF

          if(Vof % Front % cell_at_elem(c1) .ne. 0) then
            Vof % m_dot(c1) = -t_x_1 * cond_1 / 2.26e+6
          end if

          if(Vof % Front % cell_at_elem(c2) .ne. 0) then
            Vof % m_dot(c2) = -t_x_1 * cond_1 / 2.26e+6
          end if

        end if
      end do
    end if
  end do

  !-------------------!
  !   Volume source   !
  !-------------------!

  ! Here is the trick to get the sign correct:
  ! - if gas is 1 and liquid 2 => l-g =  1 => source > 0
  ! - if gas is 2 and liquid 1 => l-g = -1 => source < 0
  do c = 1, grid % n_cells
    e = Vof % Front % cell_at_elem(c)  ! Front element

    ! As Yohei and Lubomir ademantly told me - you divide with the density
    ! of the phase for "which you are solving".  And you are solving for
    ! the one which is defined as one (not zero) in the system
    if(e .ne. 0) then
      b(c) = b(c)                                                    &
           + Vof % m_dot(c) * (l-g) * Vof % Front % elem(e) % area   &
           / Vof % phase_dens(l)
    end if
  end do

  end subroutine
