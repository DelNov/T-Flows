!==============================================================================!
  subroutine Calculate_Sgs_Wale(grid)
!------------------------------------------------------------------------------!
!  Compute SGS viscosity for 'LES' by using WALE model.  
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: ONE_THIRD, TINY
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  real    :: SijdSijd(-grid % n_bnd_cells:grid % n_cells),  &
             She     (-grid % n_bnd_cells:grid % n_cells),  &
             Vor     (-grid % n_bnd_cells:grid % n_cells)
  integer :: c
  real    :: s11,  s22,  s33,  s12,  s13,  s23,  s21,  s31,  s32
  real    :: s11d, s22d, s33d, s12d, s13d, s23d, s21d, s31d, s32d
  real    :: v11,  v22,  v33,  v12,  v13,  v23,  v21,  v31,  v32
!==============================================================================!

  ! print *, '# I think there is a bug in this function (Bojan)'

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  do c = 1, grid % n_cells
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

    She(c) = 0.5 * shear(c) * shear(c)
    Vor(c) = 0.5 * vort(c) * vort(c)

    s11d =  s11*s11 + s12*s12 + s13*s13   &
         - (v11*v11 + v12*v12 + v13*v13)  &
         - ONE_THIRD * (She(c) - Vor(c))

    s22d =  s12*s12 + s22*s22 + s23*s23   &
         - (v12*v12 + v22*v22 + v23*v23)  &
         - ONE_THIRD * (She(c) - Vor(c))

    s33d =  s13*s13 + s23*s23 + s33*s33   &
         - (v13*v13 + v23*v23 + v33*v33)  &
         - ONE_THIRD * (She(c) - Vor(c))

    s12d = s11*s12 + s12*s22 + s13*s32 + (v11*v12 + v12*v22 + v13*v32) 
    s13d = s11*s13 + s12*s23 + s13*s33 + (v11*v13 + v12*v23 + v13*v33) 
    s23d = s21*s13 + s22*s23 + s23*s33 + (v21*v13 + v22*v23 + v23*v33) 

    s21d = s12d
    s31d = s13d
    s32d = s23d

    SijdSijd(c) = s11d*s11d + s22d*s22d + s33d*s33d  &
                + s12d*s12d + s13d*s13d + s23d*s23d
    
    wale_v(c) =  sqrt( abs (SijdSijd(c)**3) )          &
              / (sqrt( abs (She(c)     **5) ) +        &
                 sqrt( sqrt(SijdSijd(c)**6) ) + TINY)
  end do 

  end subroutine
