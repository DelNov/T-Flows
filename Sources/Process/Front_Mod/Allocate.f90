!==============================================================================!
  subroutine Front_Mod_Allocate(front, flow)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
  type(Field_Type), target :: flow
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  ! Take aliases to object vertex flow around
  front % pnt_flow => flow
  front % pnt_grid => flow % pnt_grid

  ! Allocate memory
  allocate(front % elem(MAX_SURFACE_ELEMENTS))
  allocate(front % vert(MAX_SURFACE_VERTICES))
  allocate(front % side(MAX_SURFACE_ELEMENTS * 3))

  ! Initialize front's local variables
  call Front_Mod_Initialize(front)

  end subroutine
