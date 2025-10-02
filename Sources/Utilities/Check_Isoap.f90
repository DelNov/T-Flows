!==============================================================================!
  program Test
!----------------------------------[Modules]-----------------------------------!
  use Polyhedron_Mod
  use Iso_Polygons_Mod
  use Isoap_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Polyhedron_Type)   :: Pol
  type(Iso_Polygons_Type) :: Iso
!-----------------------------------[Locals]-----------------------------------!
  integer       :: icellgeom, ifunc, ifile
  character(80) :: arg
!==============================================================================!

  ! Allocate Pol
  allocate(Pol % faces_n_nodes (NS));     Pol % faces_n_nodes(:)   = 0
  allocate(Pol % faces_n       (NS,NV));  Pol % faces_n      (:,:) = 0
  allocate(Pol % nodes_xyz     (NV,3));   Pol % nodes_xyz    (:,:) = 0.0
  allocate(Pol % phi           (NS));     Pol % phi          (:)   = 0.0

  ! Allocate iso-polygons
  allocate(Iso % polys_v      (NS,NV));  Iso % polys_v      (:,:) = 0
  allocate(Iso % face_index   (NV));     Iso % face_index   (:)   = 0
  allocate(Iso % polys_n_verts(NS));     Iso % polys_n_verts(:)   = 0
  allocate(Iso % verts_xyz    (NV,3));   Iso % verts_xyz    (:,:) = 0.0

  print *, '|---------------------------------------------------|'
  print *, '|      TEST PROGRAM FOR SINGLE CELLS OF ISOAP       |'
  print *, '|                                                   |'
  print *, '|             (Version 1, April 2020)               |'
  print *, '|                                                   |'
  print *, '|                       by                          |'
  print *, '|                                                   |'
  print *, '|            J. Lopez and J. Hernandez              |'
  print *, '|---------------------------------------------------|'
  print *, ''
  if(command_argument_count() .ne. 3) then
    print *, 'Correct invocation:'
    print *, ''
    print *, './test <icell_geom> <ifunc> <phi_iso>'
    print *, ''
    print *, 'Where arguments passed to test, mean the following:'
    print *, ''
    print *, 'icellgeom=  1, cube'
    print *, '         =  2, tetrahedron'
    print *, '         =  3, dodecahedron'
    print *, '         =  4, icosahedron'
    print *, '         =  5, complex cell (18 faces and 32 vertices)'
    print *, '         =101, distorted cube'
    print *, '         =102, non-convex pentagonal pyramid'
    print *, '         =103, non-convex cell obtained by subtracting a pyramid from the cube'
    print *, '         =104, stellated cube        '
    print *, '         =105, non-convex hexahedron'
    print *, '         =106, stellated dodecahedron'
    print *, '         =107, stellated icosahedron'
    print *, '         =108, hollowed cube'
    print *, '         =109, drilled cube'
    print *, '         =110, zig-zag prism'
    print *, ' ifunc   =  1, sphere with radious 0.325 centered at (0.5,0.5,0.5)'
    print *, '         =  2, torus with major radius 0.2, minor radius 0.1 and'
    print *, '               centered at (0.5,0.5,0.5)'
    print *, '         =  3, orthocircle surface centered at (1.25,1.25,1.25)'
  else
    call get_command_argument(1, arg)
    read(arg, *) icellgeom
    call get_command_argument(2, arg)
    read(arg, *) ifunc
    call get_command_argument(3, arg)
    read(arg, *) Pol % phi_iso

    call Pol % Pick_A_Test_Case(icellgeom, ifunc)

    !----------------------------------------------------------!
    !   Print the polyhedral cell to the file 'geo00000.vtk'   !
    !----------------------------------------------------------!
    ifile = 0
    call Pol % Plot_Polyhedron_Vtk("geo", ifile)

    !---------------------------------------------------------------!
    !   Iso-surface extraction: this is the core of the algorithm   !
    !---------------------------------------------------------------!
    call Isoap % Main_Isoap(Pol, Iso)

    !-------------------------------------------------------!
    !   Print the iso-polygons to the file 'iso00000.vtk'   !
    !-------------------------------------------------------!
    ifile = 0
    call Iso % Plot_Iso_Polygons_Vtk("iso", ifile)

  end if

  end program
