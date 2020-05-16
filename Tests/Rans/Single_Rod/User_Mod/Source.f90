!==============================================================================!
  subroutine User_Mod_Source(flow, phi, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Var_Type),   target :: phi
  type(Matrix_Type)        :: a_matrix
  real, dimension(:)       :: b_vector
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  type(Face_Type), pointer :: m_flux
  integer                  :: c
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  m_flux => flow % m_flux
  bulk   => flow % bulk
  u      => flow % u
  v      => flow % v
  w      => flow % w
  t      => flow % t

  !-----------------------------------------------------! 
  !                                                     !
  !   Set source depending on which variable you have   !
  !                                                     !
  !-----------------------------------------------------! 

  !-------------------------------!
  !  Set source for temperature   !
  !-------------------------------!
  if( grid % name .eq. 'FUEL') then
    if( phi % name .eq. 'T' ) then
      do c = 1, grid % n_cells
        b_vector(c) = b_vector(c) + 1.0e9 * grid % vol(c)
      end do
    end if
  end if

  end subroutine  ! fourth level comments
