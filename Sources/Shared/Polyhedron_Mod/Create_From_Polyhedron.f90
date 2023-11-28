!==============================================================================!
  subroutine Create_From_Polyhedron(Destination, Source)
!------------------------------------------------------------------------------!
!>  This subroutine is essentially as a copy constructor. It duplicates a
!>  Polyhedron_Type object (Source) into another (Destination).  This
!>  subroutine ensures that the Destination polyhedron becomes an exact
!>  copy of the Source polyhedron.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Destination  !! polyhedron being created
  type (Polyhedron_Type) :: Source       !! polyhedron being coppied
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, s
!==============================================================================!

  !----------------------------------------------------------------------!
  !   If not done yet, allocate memory for polyhedron and iso-polygons   !
  !----------------------------------------------------------------------!
  if(.not. Destination % allocated) then
    call Destination % Allocate_Polyhedron(MAX_ISOAP_FACES, MAX_ISOAP_VERTS)
  end if

  !-------------------------------------------------!
  !   (Re)initialize Destination and Iso_Polygons   !  probably not needed
  !-------------------------------------------------!
  Destination % n_nodes            = 0
  Destination % n_faces            = 0
  Destination % faces_n_nodes(:)   = 0
  Destination % faces_n      (:,:) = 0
  Destination % nodes_xyz    (:,:) = 0.0
  Destination % phi          (:)   = 0.0
  Destination % phi_iso            = 0.5
  Destination % global_node  (:)   = 0

  !-------------------------------------!
  !   Copy everything from the Source   !
  !-------------------------------------!
  Destination % n_nodes = Source % n_nodes
  Destination % n_faces = Source % n_faces
  Destination % phi_iso = Source % phi_iso
  do s = 1, Source % n_faces
    Destination % faces_n_nodes(s) = Source % faces_n_nodes(s)
    Destination % faces_n(s,:)     = Source % faces_n(s,:)
  end do
  do i = 1, Source % n_nodes
    Destination % nodes_xyz  (i,1:3) = Source % nodes_xyz  (i,1:3)
    Destination % phi        (i)     = Source % phi        (i)
    Destination % global_node(i)     = Source % global_node(i)
  end do

  end subroutine
