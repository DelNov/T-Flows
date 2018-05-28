!==============================================================================!
  subroutine Work_Mod_Allocate_Integer_Faces(grid, n)
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  integer                 :: n    ! number of arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nf
!==============================================================================!

  ! Get number of cells and boundary cells
  nf = grid % n_faces

  ! Store the pointer to the grid
  pnt_grid => grid

  ! Allocate requested memory
  allocate(i_face_01(nf));   i_face_01 = 0;   if(n .eq.  1) return
  allocate(i_face_02(nf));   i_face_02 = 0;   if(n .eq.  2) return
  allocate(i_face_03(nf));   i_face_03 = 0;   if(n .eq.  3) return
  allocate(i_face_04(nf));   i_face_04 = 0;   if(n .eq.  4) return
  allocate(i_face_05(nf));   i_face_05 = 0;   if(n .eq.  5) return
  allocate(i_face_06(nf));   i_face_06 = 0;   if(n .eq.  6) return
  allocate(i_face_07(nf));   i_face_07 = 0;   if(n .eq.  7) return
  allocate(i_face_08(nf));   i_face_08 = 0;   if(n .eq.  8) return
  allocate(i_face_09(nf));   i_face_09 = 0;   if(n .eq.  9) return
  allocate(i_face_10(nf));   i_face_10 = 0;   if(n .eq. 10) return
  allocate(i_face_11(nf));   i_face_11 = 0;   if(n .eq. 11) return
  allocate(i_face_12(nf));   i_face_12 = 0;   if(n .eq. 12) return
  allocate(i_face_13(nf));   i_face_13 = 0;   if(n .eq. 13) return
  allocate(i_face_14(nf));   i_face_14 = 0;   if(n .eq. 14) return
  allocate(i_face_15(nf));   i_face_15 = 0;   if(n .eq. 15) return
  allocate(i_face_16(nf));   i_face_16 = 0;   if(n .eq. 16) return
  allocate(i_face_17(nf));   i_face_17 = 0;   if(n .eq. 17) return
  allocate(i_face_18(nf));   i_face_18 = 0;   if(n .eq. 18) return
  allocate(i_face_19(nf));   i_face_19 = 0;   if(n .eq. 19) return
  allocate(i_face_20(nf));   i_face_20 = 0;   if(n .eq. 20) return
  allocate(i_face_21(nf));   i_face_21 = 0;   if(n .eq. 21) return
  allocate(i_face_22(nf));   i_face_22 = 0;   if(n .eq. 22) return
  allocate(i_face_23(nf));   i_face_23 = 0;   if(n .eq. 23) return
  allocate(i_face_24(nf));   i_face_24 = 0;   if(n .eq. 24) return
  allocate(i_face_25(nf));   i_face_25 = 0;   if(n .eq. 25) return
  allocate(i_face_26(nf));   i_face_26 = 0;   if(n .eq. 26) return
  allocate(i_face_27(nf));   i_face_27 = 0;   if(n .eq. 27) return
  allocate(i_face_28(nf));   i_face_28 = 0;   if(n .eq. 28) return
  allocate(i_face_29(nf));   i_face_29 = 0;   if(n .eq. 29) return
  allocate(i_face_30(nf));   i_face_30 = 0;   if(n .eq. 30) return

  end subroutine
