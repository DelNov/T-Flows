!==============================================================================!
  subroutine Control_Mod_Point_For_Monitoring_Planes(bulk, verbose)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Bulk_Mod, only: Bulk_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Bulk_Type)   :: bulk
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(3)
  real :: val(3)
!==============================================================================!

  data def / 0.0, 0.0, 0.0 /

  call Control_Mod_Read_Real_Array('POINT_FOR_MONITORING_PLANES', 3, def,  &
                                    val, verbose)
  bulk % xp = val(1)
  bulk % yp = val(2)
  bulk % zp = val(3)

  end subroutine
