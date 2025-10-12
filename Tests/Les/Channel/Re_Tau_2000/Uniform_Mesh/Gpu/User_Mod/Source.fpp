!==============================================================================!
  subroutine User_Mod_Source(Grid, Flow, Turb, phi, sc, a, b)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)          :: Grid
  type(Field_Type), target :: Flow
  type(Turb_Type)          :: Turb
  type(Var_Type),   target :: phi
  integer                  :: sc
  type(Sparse_Type)        :: a
  real, dimension(:)       :: b
!----------------------------------[Locals]------------------------------------!
  real    :: heat, area_x
  integer :: c, c2, reg, s
!==============================================================================!

  if( phi % name .eq. 'T' ) then

    !---------------------------------------------------------!
    !   If there is no bulk velocity (yet), get out of here   !
    !---------------------------------------------------------!
    if(abs(Flow % bulk % u) .lt. NANO) return

    !-------------------------------------!
    !   Compute heat flux to the domain   !
    !-------------------------------------!
    heat = 0.0
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. WALLFL) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)  ! all present
          c2 = Grid % faces_c(2,s)
          heat = heat + abs(Grid % sz(s)) * Flow % t % q(c2)
        end do
        !$tf-acc loop end

      end if
    end do

    !----------------------------------------!
    !   Cross-sectional area is hard-coded   !
    !----------------------------------------!
    area_x = 2.0 * 3.14

    !-------------------------------!
    !  Set source for temperature   !
    !-------------------------------!
    !$tf-acc loop begin
    do c = Cells_In_Domain()
      b(c) = b(c) - Flow % u % n(c) / Flow % bulk % u   &
                  * heat / area_x                       &
                  * Grid % vol(c) / Grid % per_x
    end do
    !$tf-acc loop end

  end if

  end subroutine
