!==============================================================================!
  subroutine Mass_Transfer_Estimate(Vof)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: Flow
  integer                   :: e, g, l, s, c1, c2, i_ele
  real                      :: cond_1, cond_2
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

  !----------------------------------------!
  !   Compute heat flux at the interface   !
  !----------------------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(any(Vof % Front % face_at_elem(1:2,s) .ne. 0)) then

      Vof % q_int(1,s) = 0.0
      Vof % q_int(2,s) = 0.0

      do i_ele = 1, 2
        e = Vof % Front % face_at_elem(i_ele, s)
        if(e .ne. 0) then

          ! Take conductivities from each side of the interface
          cond_1 = Vof % phase_cond(1)
          cond_2 = Vof % phase_cond(2)
          if(Vof % fun % n(c1) < 0.5) cond_1 = Vof % phase_cond(2)
          if(Vof % fun % n(c2) > 0.5) cond_2 = Vof % phase_cond(1)

          ! Compute heat fluxes from each side of the interface
          ! Units: W/(mK) * K/m * m^2 = W
          Vof % q_int(1,s) = Vof % q_int(1,s)                               &
              + cond_1 * (  Flow % t % x(c1) * Vof % Front % elem(e) % sx   &
                          + Flow % t % y(c1) * Vof % Front % elem(e) % sy   &
                          + Flow % t % z(c1) * Vof % Front % elem(e) % sz)

          Vof % q_int(2,s) = Vof % q_int(2,s)                               &
              + cond_2 * (  Flow % t % x(c2) * Vof % Front % elem(e) % sx   &
                          + Flow % t % y(c2) * Vof % Front % elem(e) % sy   &
                          + Flow % t % z(c2) * Vof % Front % elem(e) % sz)

        end if  ! e .ne. 0
      end do    ! i_ele

      ! Accumulate sources in the cells surroundig the face
      ! Unit: W * kg/J = kg / s
      if(Vof % Front % cell_at_elem(c1) .ne. 0) then
        Vof % m_dot(c1) = Vof % m_dot(c1)  &
                        - Vof % q_int(1,s) / 2.26e+6
      end if

      if(Vof % Front % cell_at_elem(c2) .ne. 0) then
        Vof % m_dot(c2) = Vof % m_dot(c2)  &
                        - Vof % q_int(1,s) / 2.26e+6
      end if

    end if      ! face is at the front

  end do

  end subroutine
