!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm,  &
                                              curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer     :: Grid
  type(Stl_Type)               :: Sphere
  integer, contiguous, pointer :: cut_count(:)
  real                         :: tot_vol, cel_vol, di
  real                         :: tet_1, tet_2, tet_3, tet_4, tet_5
  real                         :: px1, py1, pz1, px2, py2, pz2, px3, py3, pz3
  real                         :: qxi, qyi, qzi, qxj, qyj, qzj
  real                         :: xf, yf, zf, nx, ny, nz, lx, ly, lz
  integer                      :: c, s, f, i, j, i_nod, j_nod
  logical                      :: ij_cut
  integer, allocatable         :: ij_visited(:,:)
  real,    allocatable         :: x_ij_int(:)
  real,    allocatable         :: y_ij_int(:)
  real,    allocatable         :: z_ij_int(:)
!==============================================================================!

  call Work % Connect_Int_Cell(cut_count)

  ! Allocate local memory
  allocate(ij_visited(MAX_ISOAP_VERTS,MAX_ISOAP_VERTS))
  allocate(x_ij_int(MAX_ISOAP_VERTS))
  allocate(y_ij_int(MAX_ISOAP_VERTS))
  allocate(z_ij_int(MAX_ISOAP_VERTS))

  ! Take alias(es)
  Grid => Flow % pnt_grid

  !-------------------------------------------------------------!
  !   Calculate and report integral volume in traditional way   !
  !-------------------------------------------------------------!
  tot_vol = 0.0
  do c = 1, Grid % n_cells
    tot_vol = tot_vol + Grid % vol(c)
  end do
  print '(a, f12.5)', ' # Total volume from Grid is: ', tot_vol

  do c = 1, Grid % n_cells
    call Polyhedron % Extract_From_Grid(Grid, c)
    ! if(mod(c,100) .eq. 0) call Polyhedron % Plot_Polyhedron_Vtk(c)
  end do

  !---------------------------------------!
  !   Calculate individual cell volumes   !
  !---------------------------------------!
  tot_vol = 0.0
  do c = 1, Grid % n_cells
    call Polyhedron % Extract_From_Grid(Grid, c)
    call Polyhedron % Calculate_Cell_Volume(cel_vol)
    tot_vol = tot_vol + cel_vol
  end do
  print '(a, f12.5)', ' # Total volume from Polyhedron is: ', tot_vol

  !------------------------------!
  !                              !
  !   Now comes the real stuff   !
  !                              !
  !------------------------------!
  call Sphere % Create_From_File("sphere.stl")

  cut_count(:) = 0
  do c = 1, Grid % n_cells
    call Polyhedron % Extract_From_Grid(Grid, c)

    Polyhedron % phi(1:Polyhedron % n_nodes) = 0.5

    do f = 1, Sphere % n_facets

      ! STL vertex coordinates
      px1 = Sphere % x(1, f)
      py1 = Sphere % y(1, f)
      pz1 = Sphere % z(1, f)
      px2 = Sphere % x(2, f)
      py2 = Sphere % y(2, f)
      pz2 = Sphere % z(2, f)
      px3 = Sphere % x(3, f)
      py3 = Sphere % y(3, f)
      pz3 = Sphere % z(3, f)

      ! STL facet centroid
      xf = (px1 + px2 + px3) * ONE_THIRD
      yf = (py1 + py2 + py3) * ONE_THIRD
      zf = (pz1 + pz2 + pz3) * ONE_THIRD

      ! STL facet normal
      nx = Sphere % nx(f)
      ny = Sphere % ny(f)
      nz = Sphere % nz(f)

      ij_visited(:,:) = 0

      do s = 1, Polyhedron % n_faces

        ! Browse nodes in circular direction
        do i_nod = 1, Polyhedron % faces_n_nodes(s)
          j_nod = i_nod + 1; if(j_nod > Polyhedron % faces_n_nodes(s)) j_nod = 1

          i = Polyhedron % faces_n(s, i_nod)
          j = Polyhedron % faces_n(s, j_nod)

          if(ij_visited(i,j) == 0) then

          ij_cut = .false.

          qxi = Polyhedron % nodes_xyz(i, 1)
          qyi = Polyhedron % nodes_xyz(i, 2)
          qzi = Polyhedron % nodes_xyz(i, 3)

          qxj = Polyhedron % nodes_xyz(j, 1)
          qyj = Polyhedron % nodes_xyz(j, 2)
          qzj = Polyhedron % nodes_xyz(j, 3)

          ! SignedVolume(q1,p1,p2,p3)
          tet_1 = Math % Tet_Volume(qxi, qyi, qzi,   &
                                    px1, py1, pz1,   &
                                    px2, py2, pz2,   &
                                    px3, py3, pz3)

          ! SignedVolume(q2,p1,p2,p3)
          tet_2 = Math % Tet_Volume(qxj, qyj, qzj,   &
                                    px1, py1, pz1,   &
                                    px2, py2, pz2,   &
                                    px3, py3, pz3)

          ! tet_1 andn tet_2 have different signs
          if(tet_1 * tet_2 < 0.0) then

            ! SignedVolume(q1,q2,p1,p2)
            tet_3 = Math % Tet_Volume(qxi, qyi, qzi,   &
                                      qxj, qyj, qzj,   &
                                      px1, py1, pz1,   &
                                      px2, py2, pz2)

            ! SignedVolume(q1,q2,p2,p3)
            tet_4 = Math % Tet_Volume(qxi, qyi, qzi,   &
                                      qxj, qyj, qzj,   &
                                      px2, py2, pz2,   &
                                      px3, py3, pz3)
            if(tet_3 * tet_4 > 0.0) then

              ! SignedVolume(q1,q2,p3,p1)
              tet_5 = Math % Tet_Volume(qxi, qyi, qzi,   &
                                        qxj, qyj, qzj,   &
                                        px3, py3, pz3,   &
                                        px1, py1, pz1)

              if(tet_3 < 0.0 .and. tet_4 < 0.0 .and. tet_5 < 0.0) then
                ij_cut = .true.
              end if  ! tet_3, tet_4 and tet_5 have the same sign
              if(tet_3 > 0.0 .and. tet_4 > 0.0 .and. tet_5 > 0.0) then
                ij_cut = .true.
              end if  ! tet_3, tet_4 and tet_5 have the same sign

            end if  ! tet_3 and tet_4 have the same sign

          end if  ! tet_1 and tet_2 have different signs

          ! There was a cut_count in between i and j
          if(ij_cut) then
            cut_count(c) = cut_count(c) + 1

            ! Vector connecting nodes i and j
            lx = qxj - qxi
            ly = qyj - qyi
            lz = qzj - qzi

            ! Distance from i to a point on the plane
            ! (I picked center of the facet)
            di = (  (xf - qxi) * nx     &
                  + (yf - qyi) * ny     &
                  + (zf - qzi) * nz  )  &
               / (lx * nx + ly * ny + lz * nz)

            ! Intersection point
            x_ij_int(cut_count(c)) = qxi + di * lx
            y_ij_int(cut_count(c)) = qyi + di * ly
            z_ij_int(cut_count(c)) = qzi + di * lz

            ! Mark nodes which are on either side of interface ...
            ! ... using the scalar product with STL facet centre
            if( nx*(qxi-xf) + ny*(qyi-yf) + nz*(qzi-zf) > 0.0 ) then
              Polyhedron % phi(i) = 1.0
            else
              Polyhedron % phi(i) = 0.0
            end if
            if( nx*(qxj-xf) + ny*(qyj-yf) + nz*(qzj-zf) > 0.0 ) then
              Polyhedron % phi(j) = 1.0
            else
              Polyhedron % phi(j) = 0.0
            end if

            ij_visited(i,j) = cut_count(c)
            ij_visited(j,i) = cut_count(c)
          end if   ! .not. visited

          end if  ! ij_cut

        end do  ! through nodes of the face

      end do    ! through polyhedron faces

    end do  ! through STL facets

    if(cut_count(c) .ge. 1) call Polyhedron % Plot_Polyhedron_Vtk(c)
  end do  ! through cells

  call Work % Disconnect_Int_Cell(cut_count)

  ! Exit gracefully
  call Comm_Mod_End
  stop

  end subroutine
