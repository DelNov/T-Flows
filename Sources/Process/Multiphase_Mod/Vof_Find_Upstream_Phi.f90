!==============================================================================!
subroutine Multiphase_Mod_Vof_Find_Upstream_phi(phi, phi_x, phi_y, phi_z,  &
                                                s, donor, accept, phi_u)
!------------------------------------------------------------------------------!
!   Computes the value of phi at a imaginary upstream cell. This is based on   !
!   Work of Zhang (2014) "Assessment of different reconstruction techniques    !
!          for implementing the NVSF schemes on unstructured meshes"           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)           :: phi
  real                     :: phi_x(-phi % pnt_grid % n_bnd_cells:  &
                                     phi % pnt_grid % n_cells),     &
                              phi_y(-phi % pnt_grid % n_bnd_cells:  &
                                     phi % pnt_grid % n_cells),     &
                              phi_z(-phi % pnt_grid % n_bnd_cells:  &
                                     phi % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real,            pointer :: sx(:), sy(:), sz(:), dx(:), dy(:), dz(:)
  real                     :: du_orig, du_pred, dotprod, du, dd
  real                     :: dd_vect(3), du_vect(3), signo, phi_u
  real                     :: min_phi, max_phi
  integer                  :: s, ss, idonor, n_nodes_c, n_node_f, n1, n2, nod
  integer                  :: c, c1, c2, donor, accept, c1_glo, c2_glo
  logical                  :: out_face
!==============================================================================!

  ! Take aliases
  grid => phi % pnt_grid
  sx   => grid % sx;  sy   => grid % sy;  sz   => grid % sz
  dx   => grid % dx;  dy   => grid % dy;  dz   => grid % dz

  !find min and max phi:

  min_phi =  HUGE
  max_phi = -HUGE
  do idonor = 1, grid % cells_n_faces(donor)

    ss = grid % cells_f(idonor, donor)
    c1 = grid % faces_c(1,ss)
    c2 = grid % faces_c(2,ss)

    if (c1 .ne. accept) then
      min_phi = min(min_phi, phi % n(c1))
      max_phi = max(max_phi, phi % n(c1))
    end if

    if (c2 .ne. accept) then
      min_phi = min(min_phi, phi % n(c2))
      max_phi = max(max_phi, phi % n(c2))
    end if
  end do

  if (donor == grid % faces_c(1,s)) then
    signo = 1.0
  else
    signo = -1.0
  end if

  du_orig = norm2((/dx(s),dy(s),dz(s)/))

  dd = abs(dot_product((/dx(s),dy(s),dz(s)/),                       &
                       (/sx(s),sy(s),sz(s)/) / grid % s(s)))

  dd_vect = dd * (/sx(s),sy(s),sz(s)/) / grid % s(s) * signo

  du_pred = 0.0

  n_nodes_c = grid % cells_n_nodes(donor)

  do n1 = 1, grid % cells_n_nodes(donor)
    out_face = .true.
    loop_face: do n2 = 1, grid % faces_n_nodes(s)
      if (grid % cells_n(n1,donor) == grid % faces_n(n2,s)) then
        n_nodes_c = n_nodes_c - 1 
        out_face = .false.
        exit loop_face
      end if
    end do loop_face

    if (out_face) then  ! node doesn't belong to face s
      dotprod = dot_product(                                                 &
                (/grid % xn(grid % cells_n(n1,donor)) - grid % xc(donor),    &
                  grid % yn(grid % cells_n(n1,donor)) - grid % yc(donor),    &
                  grid % zn(grid % cells_n(n1,donor)) - grid % zc(donor)/),  &
                  (/sx(s),sy(s),sz(s)/) / grid % s(s) * signo) 
      if (dotprod > 0.0) then
        du_pred = du_pred + dotprod
      else
        n_nodes_c = n_nodes_c - 1
      end if
    end if
  end do

  du_pred = du_pred / real(n_nodes_c)

  du = min(du_pred, du_orig)

  du_vect = - du * (/sx(s),sy(s),sz(s)/) / grid % s(s) * signo

  phi_u = (phi % n(accept) - phi % n(donor) - ( phi_x(donor) * dd_vect(1)     &
                                              + phi_y(donor) * dd_vect(2)     &
                                              + phi_z(donor) * dd_vect(3) ))  &
        * du * du / (dd * dd) + ( phi_x(donor) * du_vect(1)                   &
                                + phi_y(donor) * du_vect(2)                   &
                                + phi_z(donor) * du_vect(3) )                 &
        + phi % n(donor)

  phi_u = min(max(phi_u, min_phi), max_phi)

  end subroutine
