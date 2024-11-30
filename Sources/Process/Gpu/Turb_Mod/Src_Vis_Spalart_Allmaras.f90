!==============================================================================!
  subroutine Src_Vis_Spalart_Allmaras(Turb, Grid, Flow)
!------------------------------------------------------------------------------!
!   Computes the source terms in vis transport equation.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),      target :: Turb
  type(Grid_Type)               :: Grid
  type(Field_Type),      target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  real,    contiguous, pointer :: vis_x(:), vis_y(:), vis_z(:)
  real,    contiguous, pointer :: val(:)
  integer, contiguous, pointer :: dia(:), pos(:,:)
  real,    contiguous, pointer :: b(:)
  integer                      :: c
  real                         :: x_rat, f_v1, f_v2, f_w, ss
  real                         :: dist_v, prod_v, r, gg, dif, dist
!==============================================================================!

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % C % dia
  pos => Flow % Nat % C % pos
  b   => Flow % Nat % b

  call Grad_Variable(Flow, Grid, Turb % vis)

  vis_x => Flow % phi_x
  vis_y => Flow % phi_y
  vis_z => Flow % phi_z

# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("vis_gradients",        &
                               vector_name = "grads",  &
                               vector_cell = (/vis_x, vis_y, vis_z/))
# endif

  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then

    do c = Cells_In_Domain_And_Buffers()

      dist = Grid % wall_dist(c)
      if(Turb % model .eq. DES_SPALART) then
        ! What is 0.65 here?  A ghost number
        dist = min(Grid % wall_dist(c), 0.65 * Turb % h_max(c))
      end if

      !---------------------------------!
      !   Compute the production term   !
      !---------------------------------!
      x_rat  = Turb % vis % n(c) / (Flow % viscosity(c) / Flow % density(c))
      f_v1   = x_rat**3 / (x_rat**3 + Turb % c_v1**3)
      f_v2   = 1.0 - x_rat/(1.0 + x_rat*f_v1)
      ss     = Flow % vort(c)   &
             + Turb % vis % n(c) * f_v2 / (Turb % kappa**2 * dist**2)
      prod_v = Turb % c_b1 * Flow % density(c) * ss * Turb % vis % n(c)
      b(c)   = b(c) + prod_v * Grid % vol(c)

      !----------------------------------!
      !   Compute the destruction term   !
      !----------------------------------!
      r      = Turb % vis % n(c) / (ss * Turb % kappa**2 * dist**2)
      gg     = r + Turb % c_w2*(r**6 - r)
      f_w    = gg*((1.0 + Turb % c_w3**6)  &
             / (gg**6 + Turb % c_w3**6))**ONE_SIXTH
      dist_v = Turb % c_w1 * Flow % density(c) * f_w  &
             * (Turb % vis % n(c) / dist**2)
      val(dia(c)) = val(dia(c)) + dist_v * Grid % vol(c)

      !--------------------------------------------!
      !   Compute the first-order diffusion term   !
      !--------------------------------------------!
      dif   = Turb % c_b2                                &
            * Flow % density(c)                          &
            * (vis_x(c)**2 + vis_y(c)**2 + vis_z(c)**2)  &
            / Turb % vis % sigma
      b(c)  = b(c) + dif * Grid % vol(c)

    end do

  end if

  end subroutine
