!==============================================================================!
  subroutine Main_Vof(Vof, Flow, turb, Sol, curr_dt)
!------------------------------------------------------------------------------!
!   Initialize inteface tracking simulations                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),     target :: Vof
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Solver_Type)           :: Sol
  integer, intent(in)         :: curr_dt     ! time step
!==============================================================================!

  if(Flow % with_interface) then

    !-----------------!
    !                 !
    !   Compute vof   !
    !                 !
    !-----------------!

    !--------------------------------!
    !   Advance vof function (fun)   !
    !--------------------------------!
    call Update_Boundary_Values(Flow, turb, Vof, 'VOF')
    call Vof % Compute_Vof(Sol, Flow % dt, curr_dt)

    !------------------------------------------------!
    !   Prepare smooth variant of the vof function   !
    !    for computation of normals and curvature    !
    !------------------------------------------------!
    if(Vof % surface_tension > TINY) then
      call Vof % Smooth_For_Curvature_Csf()
      call Vof % Curvature_Csf()
    end if

    !---------------------------------------------------------------!
    !   Update properties for other conservation equations to use   !
    !    (Maybe redundant since it is called from Main_Pro too.)    !
    !---------------------------------------------------------------!
    call Vof % Update_Physical_Properties()

    !-------------------------------!
    !                               !
    !   Create front or a surface   !
    !                               !
    !-------------------------------!

    !---------------------!
    !   Created a front   !
    !---------------------!
    if(Vof % track_front) then
      call Vof % Front % Place_Front_At_Value(Vof % fun,     &
                                              Vof % smooth,  &
                                              0.5,           &
                                              .true.)  ! don't print messages
      call Vof % Front % Print_Front_Statistics()
    end if

    !----------------------!
    !   Create a surface   !
    !----------------------!
    if(Vof % track_surface) then
      call Vof % Surf % Place_Surf_At_Value(Vof % fun,     &
                                            Vof % smooth,  &
                                            0.5,           &
                                            .true.)  ! don't print messages
      call Vof % Surf % Improve_Mesh_Quality(Vof % smooth,  &
                                             0.5,           &
                                            .true.)  ! don't print messages
!     call Vof % Surf % Calculate_Curvatures_From_Elems()
      call Vof % Surf % Print_Front_Statistics()
    end if

  end if

  end subroutine
