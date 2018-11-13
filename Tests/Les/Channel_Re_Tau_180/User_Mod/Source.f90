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
  ! |dT/dx * Ux * volume 
  ! |
  !/
  ! dT/dx is derived from condition that there is no
  ! change of energy in the system. It means
  ! mass_flux * cp * dT - Q * Area = 0, Area = B*dx*Nwall,
  !  dT/dx = Q*B*Nwall/(mass_flux*cp),
  ! where Q is heat flux through the wall, B is
  ! channel width and Nwall is number of heated walls.  

  !-------------------------------!
  !  Set source for temperature   !
  !-------------------------------!
  if( phi % name .eq. 'T' ) then  
    do c = 1, grid % n_cells
      
      b_vector(c) = b_vector(c)   &
                  - 3.14 * 2* heat_flux * u % n(c) &
      / (bulk % flux_x * capacity) * grid % vol(c)
    end do
  end if

  end subroutine  ! fourth level comments
