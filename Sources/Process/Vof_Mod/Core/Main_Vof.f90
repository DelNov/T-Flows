!==============================================================================!
  subroutine Main_Vof(Vof, Flow, Turb, Sol)
!------------------------------------------------------------------------------!
!   Initialize inteface tracking simulations                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),     target :: Vof
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Solver_Type)           :: Sol
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  if(Flow % with_interface) then

    !-----------------!
    !                 !
    !   Compute vof   !
    !                 !
    !-----------------!

    !------------------------------------!
    !   Advance vof function (fun) and   !
    !    re-create its smooth variant    !
    !------------------------------------!
    call Vof % Compute_Vof(Sol, Flow % dt)

    !------------------------------------------------!
    !   Prepare smooth variant of the vof function   !
    !    for computation of normals and curvature    !
    !------------------------------------------------!
    call Vof % Smooth_Vof_And_Compute_Surface_Normals()
    call Vof % Curvature_Csf()

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

    !--------------------!
    !   Create a front   !
    !--------------------!
    if(Vof % track_front) then
      call Vof % Front % Place_Front_At_Value(Vof % fun,     &
                                              .false.)  ! don't print messages
      call Vof % Front % Print_Front_Statistics()
    end if

    !----------------------!
    !   Create a surface   !
    !----------------------!
    if(Vof % track_surface) then
      call Vof % Surf % Place_Surf_At_Value(Vof % fun,     &
                                            Vof % smooth,  &
                                            .true.)  ! don't print messages
      call Vof % Surf % Improve_Mesh_Quality(Vof % smooth)
!     call Vof % Surf % Calculate_Curvatures_From_Elems()
      call Vof % Surf % Print_Front_Statistics()
    end if

  end if

  end subroutine
