!==============================================================================!
  subroutine Work_Mod_Allocate_Integer_Nodes(grid, n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid(:)
  integer         :: n    ! number of arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nn
!==============================================================================!

  if(n .eq. 0) return

  ! Get number of nodes
  nn = maxval(grid(1:size(grid)) % n_nodes)

  ! Allocate requested memory
  allocate(i_node_01(nn));  i_node_01 = 0;  if(n .eq.  1) return
  allocate(i_node_02(nn));  i_node_02 = 0;  if(n .eq.  2) return
  allocate(i_node_03(nn));  i_node_03 = 0;  if(n .eq.  3) return
  allocate(i_node_04(nn));  i_node_04 = 0;  if(n .eq.  4) return
  allocate(i_node_05(nn));  i_node_05 = 0;  if(n .eq.  5) return
  allocate(i_node_06(nn));  i_node_06 = 0;  if(n .eq.  6) return
  allocate(i_node_07(nn));  i_node_07 = 0;  if(n .eq.  7) return
  allocate(i_node_08(nn));  i_node_08 = 0;  if(n .eq.  8) return
  allocate(i_node_09(nn));  i_node_09 = 0;  if(n .eq.  9) return
  allocate(i_node_10(nn));  i_node_10 = 0;  if(n .eq. 10) return
  allocate(i_node_11(nn));  i_node_11 = 0;  if(n .eq. 11) return
  allocate(i_node_12(nn));  i_node_12 = 0;  if(n .eq. 12) return

  end subroutine
