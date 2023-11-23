!==============================================================================!
  subroutine User_Mod_Force(Flow, Por, ui, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for velocity.      !
!   It is called from "Compute_Velocity" function, just before calling the     !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)    :: Flow
  type(Porosity_Type) :: Por
  type(Var_Type)      :: ui        ! velocity component
  type(Matrix_Type)   :: a_matrix  ! system matrix
  real, dimension(:)  :: b_vector  ! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c, reg
  real                      :: u_mag, a_coef
!==============================================================================!

  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  !----------------------------------------------------!
  !   Set source depending on the velocity component   !
  !----------------------------------------------------!
  do reg = 1, Por % n_regions
    do c = 1, Grid % n_cells
      if(Grid % por(c) .eq. reg) then

        u_mag = sqrt(u % n(c)**2 + v % n(c)**2 + w % n(c)**2)

        if(Grid % wall_dist(c) < 20.0) then

          !-------------------------------------------------------!
          !   The a_coef is approximated to follow the coeff.     !
          !   used by Tom Grylls and Maarten van Reeuwijk in      !
          !   publication "Tree model with drag, transpiration,   !
          !   shading and deposition: Identification of cooling   !
          !   regimes and large-eddy simulation"                  !
          !-------------------------------------------------------!
          a_coef = 1.65/8.172                                                &
                   * (3.5642635513332240e-006 * Grid % wall_dist(c)**5       &
                    - 1.3555088944116367e-004 * Grid % wall_dist(c)**4       &
                    + 1.7125968146226682e-004 * Grid % wall_dist(c)**3       &
                    + 2.0848654211247653e-002 * Grid % wall_dist(c)**2       &
                    + 4.2234984220098151e-003 * Grid % wall_dist(c)          &
                    + 5.0417568827710268e-001)
        end if

        if( ui % name .eq. 'U' ) then
          b_vector(c) = b_vector(c)  &
                      - 0.5 * 0.15 * a_coef * u_mag * u % n(c) * Grid % vol(c)
        else if( ui % name .eq. 'V' ) then
          b_vector(c) = b_vector(c)  &
                      - 0.5 * 0.15 * a_coef * u_mag * v % n(c) * Grid % vol(c)
        else if( ui % name .eq. 'W' ) then
          b_vector(c) = b_vector(c)  &
                      - 0.5 * 0.15 * a_coef * u_mag * w % n(c) * Grid % vol(c)
        end if
      end if
    end do
  end do

  end subroutine
