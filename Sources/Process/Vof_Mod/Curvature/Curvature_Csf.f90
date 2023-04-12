!==============================================================================!
  subroutine Curvature_Csf(Vof)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF using Gauss theorem        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: fun
  type(Var_Type),   pointer :: smooth
  integer                   :: c, c1, c2, s, nb, nc
  real                      :: v1(3), v2(3), v3(3), v4(3)
  real                      :: norm_grad, dotprod, mag_ntan, theta
  real                      :: nwallx, nwally, nwallz
  real                      :: nintx, ninty, nintz, ntanx, ntany, ntanz
  real, contiguous, pointer :: div_x(:), div_y(:), div_z(:)
!==============================================================================!

  call Work % Connect_Real_Cell(div_x, div_y, div_z)

  Grid   => Vof % pnt_grid
  Flow   => Vof % pnt_flow
  fun    => Vof % fun
  smooth => Vof % smooth

  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  Vof % curv = 0.0

  !-------------------------------!
  !   Normalize vector at cells   !
  !-------------------------------!
  do c = 1, Grid % n_cells
    norm_grad = sqrt(  smooth % x(c) ** 2  &
                     + smooth % y(c) ** 2  &
                     + smooth % z(c) ** 2)
    if(norm_grad >= FEMTO) then
      Vof % nx(c) = smooth % x(c) / norm_grad
      Vof % ny(c) = smooth % y(c) / norm_grad
      Vof % nz(c) = smooth % z(c) / norm_grad
    else
      Vof % nx(c) = 0.0
      Vof % ny(c) = 0.0
      Vof % nz(c) = 0.0
    end if
  end do

  call Grid % Exchange_Cells_Real(Vof % nx(-nb:nc))
  call Grid % Exchange_Cells_Real(Vof % ny(-nb:nc))
  call Grid % Exchange_Cells_Real(Vof % nz(-nb:nc))

  !---------------------------------------!
  !   Tangent vector to walls/symmetries  !
  !---------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL   .or.   &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL .or.   &
         Grid % Bnd_Cond_Type(c2) .eq. SYMMETRY) then

        norm_grad = sqrt(  Vof % nx(c1)**2   &
                         + Vof % ny(c1)**2   &
                         + Vof % nz(c1)**2)

        if(norm_grad > FEMTO) then
          v1(1) = Grid % sx(s)
          v1(2) = Grid % sy(s)
          v1(3) = Grid % sz(s)
          v2(1) = Vof % nx(c1)
          v2(2) = Vof % ny(c1)
          v2(3) = Vof % nz(c1)
          v3 = Math % Cross_Product(v1, v2)
          v4 = Math % Cross_Product(v3, v1)

          ! Projection on v4
          norm_grad = sqrt(v4(1)**2 + v4(2)**2 + v4(3)**2)

          if(norm_grad > FEMTO) then
            Vof % nx(c2) = v4(1) / norm_grad
            Vof % ny(c2) = v4(2) / norm_grad
            Vof % nz(c2) = v4(3) / norm_grad
          end if
        end if
      else
        Vof % nx(c2) = Vof % nx(c1)
        Vof % ny(c2) = Vof % ny(c1)
        Vof % nz(c2) = Vof % nz(c1)
      end if

    end if  ! c2 < 0
  end do

  !----------------------------------------!
  !   Correct for contact angle at walls   !
  !----------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.    &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        ! Accumulate values of faces
        norm_grad = sqrt(  Vof % nx(c1)**2   &
                         + Vof % ny(c1)**2   &
                         + Vof % nz(c1)**2)

        if(norm_grad > FEMTO) then

          ! unit normal vector of interface
          nintx = Vof % nx(c1)
          ninty = Vof % ny(c1)
          nintz = Vof % nz(c1)

          ! unit wall normal vector directed into wall
          nwallx = Grid % sx(s)/Grid %s(s)
          nwally = Grid % sy(s)/Grid %s(s)
          nwallz = Grid % sz(s)/Grid %s(s)

          ! ntan vector lies in wall and is normal to contact line, see Brackbill's paper
          ! ntan = (nint - (nint.nwall) nwall) / abs(nint - (nint.nwall) nwall)
          dotprod = nintx * nwallx + ninty * nwally + nintz * nwallz
          ntanx = nintx - dotprod * nwallx
          ntany = ninty - dotprod * nwally
          ntanz = nintz - dotprod * nwallz
          ! normalize
          mag_ntan = sqrt( ntanx**2.0 + ntany**2.0 + ntanz**2.0 + FEMTO)
          ntanx = ntanx/mag_ntan
          ntany = ntany/mag_ntan
          ntanz = ntanz/mag_ntan

          theta = fun % q(c2) * PI /180.0
          Vof % nx(c2) = nwallx * cos(theta) + ntanx * sin(theta)
          Vof % ny(c2) = nwally * cos(theta) + ntany * sin(theta)
          Vof % nz(c2) = nwallz * cos(theta) + ntanz * sin(theta)

          ! debug
          !IF (abs(Grid % zf(s)-0.0e-3).LT.1e-6) THEN
          !  WRITE(*,*)'Curv:x',Grid % xf(s), Grid % yf(s), Grid % zf(s)
          !  WRITE(*,*)'Curv:d',Grid % dx(s), Grid % dy(s), Grid % dz(s)
          !  WRITE(*,*)'Curv:s',Grid % sx(s), Grid % sy(s), Grid % sz(s), Grid % s(s) 
          !ENDIF

        end if
      end if  ! if WALL
    end if  ! c2 < 0

  end do

  call Grid % Exchange_Cells_Real(Vof % nx(-nb:nc))
  call Grid % Exchange_Cells_Real(Vof % ny(-nb:nc))
  call Grid % Exchange_Cells_Real(Vof % nz(-nb:nc))

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  ! Find divergence of normals
  call Flow % Grad_Component(Vof % nx(-nb:nc), 1, div_x(-nb:nc))
  call Flow % Grad_Component(Vof % ny(-nb:nc), 2, div_y(-nb:nc))
  call Flow % Grad_Component(Vof % nz(-nb:nc), 3, div_z(-nb:nc))

  Vof % curv(-nb:nc) = Vof % curv(-nb:nc) - div_x(-nb:nc)
  Vof % curv(-nb:nc) = Vof % curv(-nb:nc) - div_y(-nb:nc)
  Vof % curv(-nb:nc) = Vof % curv(-nb:nc) - div_z(-nb:nc)

  call Grid % Exchange_Cells_Real(Vof % curv)

  ! call Grid % Save_Debug_Vtu('curv_sharp',scalar_cell=Vof % curv,  &
  !                                         scalar_name='curv_sharp')

  call Vof % Smooth_Curvature()

  ! call Grid% Save_Debug_Vtu('curv_smooth',scalar_cell=Vof % curv,  &
  !                                         scalar_name='curv_smooth')

  do c = 1, Grid % n_cells
    if(smooth % n(c) < 0.01 .or. smooth % n(c) > 0.99) Vof % curv(c) = 0.0
  end do

  call Work % Disconnect_Real_Cell(div_x, div_y, div_z)

  end subroutine
