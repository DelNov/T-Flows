!==============================================================================!
  subroutine User_Mod_Save_Results(flow, turb, save_name)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  character(len=*)         :: save_name
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

  end subroutine
