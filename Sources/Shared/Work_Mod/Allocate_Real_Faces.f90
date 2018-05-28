!==============================================================================!
  subroutine Work_Mod_Allocate_Real_Faces(grid, n)
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
  allocate(r_face_01(nf));   r_face_01 = 0.0;   if(n .eq.  1) return
  allocate(r_face_02(nf));   r_face_02 = 0.0;   if(n .eq.  2) return
  allocate(r_face_03(nf));   r_face_03 = 0.0;   if(n .eq.  3) return
  allocate(r_face_04(nf));   r_face_04 = 0.0;   if(n .eq.  4) return
  allocate(r_face_05(nf));   r_face_05 = 0.0;   if(n .eq.  5) return
  allocate(r_face_06(nf));   r_face_06 = 0.0;   if(n .eq.  6) return
  allocate(r_face_07(nf));   r_face_07 = 0.0;   if(n .eq.  7) return
  allocate(r_face_08(nf));   r_face_08 = 0.0;   if(n .eq.  8) return
  allocate(r_face_09(nf));   r_face_09 = 0.0;   if(n .eq.  9) return
  allocate(r_face_10(nf));   r_face_10 = 0.0;   if(n .eq. 10) return
  allocate(r_face_11(nf));   r_face_11 = 0.0;   if(n .eq. 11) return
  allocate(r_face_12(nf));   r_face_12 = 0.0;   if(n .eq. 12) return
  allocate(r_face_13(nf));   r_face_13 = 0.0;   if(n .eq. 13) return
  allocate(r_face_14(nf));   r_face_14 = 0.0;   if(n .eq. 14) return
  allocate(r_face_15(nf));   r_face_15 = 0.0;   if(n .eq. 15) return
  allocate(r_face_16(nf));   r_face_16 = 0.0;   if(n .eq. 16) return
  allocate(r_face_17(nf));   r_face_17 = 0.0;   if(n .eq. 17) return
  allocate(r_face_18(nf));   r_face_18 = 0.0;   if(n .eq. 18) return
  allocate(r_face_19(nf));   r_face_19 = 0.0;   if(n .eq. 19) return
  allocate(r_face_20(nf));   r_face_20 = 0.0;   if(n .eq. 20) return
  allocate(r_face_21(nf));   r_face_21 = 0.0;   if(n .eq. 21) return
  allocate(r_face_22(nf));   r_face_22 = 0.0;   if(n .eq. 22) return
  allocate(r_face_23(nf));   r_face_23 = 0.0;   if(n .eq. 23) return
  allocate(r_face_24(nf));   r_face_24 = 0.0;   if(n .eq. 24) return
  allocate(r_face_25(nf));   r_face_25 = 0.0;   if(n .eq. 25) return
  allocate(r_face_26(nf));   r_face_26 = 0.0;   if(n .eq. 26) return
  allocate(r_face_27(nf));   r_face_27 = 0.0;   if(n .eq. 27) return
  allocate(r_face_28(nf));   r_face_28 = 0.0;   if(n .eq. 28) return
  allocate(r_face_29(nf));   r_face_29 = 0.0;   if(n .eq. 29) return
  allocate(r_face_30(nf));   r_face_30 = 0.0;   if(n .eq. 30) return

  end subroutine
