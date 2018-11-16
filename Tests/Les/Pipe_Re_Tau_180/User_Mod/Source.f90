!==============================================================================!
  subroutine User_Mod_Source(grid, phi, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Var_Mod
  use Matrix_Mod
  use Flow_Mod
  use Const_Mod, only: PI
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Locals]------------------------------------!
  integer :: c
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)    :: grid
  type(Var_Type)     :: phi
  type(Matrix_Type)  :: a_matrix
  real, dimension(:) :: b_vector
!==============================================================================!

  !-----------------------------------------------------! 
  !                                                     !
  !   Set source depending on which variable you have   !
  !                                                     !
  !-----------------------------------------------------! 
  ! 
  !
  !  /
  ! |
  ! |dT/dz * Uz * volume 
  ! |
  !/
  ! dT/dz is derived from condition that there is no
  ! change of energy in the system. It means
  ! mass_flux * cp * dT - Q * Area = 0, Area = 2 r pi * dz ,
  !  dT/dz = 2 r pi Q / (mass_flux*cp). 
  ! where Q is heat flux through the wall.

  !-------------------------------!
  !  Set source for temperature   !
  !-------------------------------!
  if( phi % name .eq. 'T' ) then  
    do c = 1, grid % n_cells
      
      b_vector(c) = b_vector(c)   &
                  - 2.0 * pi * heat_flux * w % n(c) &
      / (bulk % flux_z * capacity) * grid % vol(c)

    end do
  end if

  end subroutine  ! fourth level comments
