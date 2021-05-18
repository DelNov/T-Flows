!==============================================================================!
  subroutine Work_Mod_Allocate_Real_Nodes(Grid, n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid(:)
  integer         :: n    ! number of arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nn
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of nodes
  nn = maxval(Grid(1:size(Grid)) % n_nodes)

  ! Allocate requested memory
  allocate(r_node_01(nn));  r_node_01 = 0.0;  if(n .eq.  1) return
  allocate(r_node_02(nn));  r_node_02 = 0.0;  if(n .eq.  2) return
  allocate(r_node_03(nn));  r_node_03 = 0.0;  if(n .eq.  3) return
  allocate(r_node_04(nn));  r_node_04 = 0.0;  if(n .eq.  4) return
  allocate(r_node_05(nn));  r_node_05 = 0.0;  if(n .eq.  5) return
  allocate(r_node_06(nn));  r_node_06 = 0.0;  if(n .eq.  6) return
  allocate(r_node_07(nn));  r_node_07 = 0.0;  if(n .eq.  7) return
  allocate(r_node_08(nn));  r_node_08 = 0.0;  if(n .eq.  8) return
  allocate(r_node_09(nn));  r_node_09 = 0.0;  if(n .eq.  9) return
  allocate(r_node_10(nn));  r_node_10 = 0.0;  if(n .eq. 10) return
  allocate(r_node_11(nn));  r_node_11 = 0.0;  if(n .eq. 11) return
  allocate(r_node_12(nn));  r_node_12 = 0.0;  if(n .eq. 12) return

  end subroutine
