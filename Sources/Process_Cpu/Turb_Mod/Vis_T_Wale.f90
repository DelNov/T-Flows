!==============================================================================!
  subroutine Vis_T_Wale(Turb)
!------------------------------------------------------------------------------!
!  Compute SGS viscosity for 'LES' by using LES_WALE model.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c
  real                      :: s11, s22, s33,  s12, s13, s23,  s21, s31, s32
  real                      :: s11d,s22d,s33d, s12d,s13d,s23d, s21d,s31d,s32d
  real                      :: v11, v22, v33,  v12, v13, v23,  v21, v31, v32
  real, contiguous, pointer :: sijd_sijd(:), shear2(:), vort2(:)
!==============================================================================!

  call Work % Connect_Real_Cell(sijd_sijd, shear2, vort2)

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  do c = Cells_In_Domain_And_Buffers()
    s11 = u % x(c)
    s22 = v % y(c)
    s33 = w % z(c)
    s12 = 0.5*(v % x(c) + u % y(c))
    s13 = 0.5*(u % z(c) + w % x(c))
    s23 = 0.5*(v % y(c) + w % y(c))
    s21 = s12
    s31 = s13
    s32 = s23

    v11 = 0.0
    v22 = 0.0
    v33 = 0.0
    v12 = 0.5*(v % x(c) - u % y(c))
    v13 = 0.5*(u % z(c) - w % x(c))
    v23 = 0.5*(v % y(c) - w % y(c))
    v21 = -v12
    v31 = -v13
    v32 = -v23

    shear2(c) = 0.5 * Flow % shear(c) * Flow % shear(c)
    vort2(c)  = 0.5 * Flow % vort(c) * Flow % vort(c)

    s11d =  s11*s11 + s12*s12 + s13*s13   &
         - (v11*v11 + v12*v12 + v13*v13)  &
         - ONE_THIRD * (shear2(c) - vort2(c))

    s22d =  s12*s12 + s22*s22 + s23*s23   &
         - (v12*v12 + v22*v22 + v23*v23)  &
         - ONE_THIRD * (shear2(c) - vort2(c))

    s33d =  s13*s13 + s23*s23 + s33*s33   &
         - (v13*v13 + v23*v23 + v33*v33)  &
         - ONE_THIRD * (shear2(c) - vort2(c))

    s12d = s11*s12 + s12*s22 + s13*s32 + (v11*v12 + v12*v22 + v13*v32) 
    s13d = s11*s13 + s12*s23 + s13*s33 + (v11*v13 + v12*v23 + v13*v33) 
    s23d = s21*s13 + s22*s23 + s23*s33 + (v21*v13 + v22*v23 + v23*v33) 

    s21d = s12d
    s31d = s13d
    s32d = s23d

    sijd_sijd(c) = s11d*s11d + s22d*s22d + s33d*s33d  &
                 + s12d*s12d + s13d*s13d + s23d*s23d

    Turb % wale_v(c) =  sqrt( abs (sijd_sijd(c)**3) )   &
              / (sqrt( abs (shear2(c)   **5) ) +        &
                 sqrt( sqrt(sijd_sijd(c)**6) ) + TINY)
  end do

  call Work % Disconnect_Real_Cell(sijd_sijd, shear2, vort2)

  end subroutine
