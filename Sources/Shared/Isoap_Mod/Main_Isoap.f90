!---------------------------------------------------------------------!
!                               ISOAP                                 !
!... Isosurface extraction on a general polyhedron                    !
!---------------------------------------------------------------------!
!>  Original main Isoap function. Calls Isopol to arrange iso-vertices.
!>  No modifications have been made to the original code from the
!>  authors to align with T-Flows coding conventions. This decision is
!>  to maintain compatibility with potential future updates or
!>  corrections from the original authors. Given the complex nature of
!>  these algorithms, which are not entirely within our current domain
!>  of expertise, preserving the original structure and logic is deemed
!>  the most prudent approach.
!---------------------------------------------------------------------!
! On entry:                                                           !
!==========                                                           !
! IPV      = array containing the indices of the vertices of each     !
!            face of the polyhedron                                   !
! NIPV     = number of vertices of each face of the polyhedron        !
! NTP      = number of vertices of the polyhedron                     !
! NTS      = number of faces of the polyhedron                        !
! PHI      = scalar value assigned to each vertex of the polyhedron   !      
! PHIISO   = iso value                                                !
! VERTP    = coordinates of the vertices of the polyhedron            !
! On return:                                                          !
!===========                                                          !
! IPVISO   = array conatining the indices of the iso-vertices of each !
!            iso-polygon                                              !
! ISOEFACE = face index of the polyhedron over which is constructed   !
!            each iso-edge                                            !
! NIPVISO  = number of iso-vertices of each iso-polygon               !
! NISO     = number of iso-polygons                                   !
! VERTISO  = coordinates of the iso vertices                          !      
!---------------------------------------------------------------------!
  SUBROUTINE MAIN_ISOAP(ISOAP, POLYHEDRON, ISO_POLYGONS)
!---------------------------------------------------------------------!
  IMPLICIT NONE
!---------------------------------------------------------------------!
  CLASS(ISOAP_TYPE)               :: ISOAP
  TYPE(POLYHEDRON_TYPE),   TARGET :: POLYHEDRON
  TYPE(ISO_POLYGONS_TYPE), TARGET :: ISO_POLYGONS
!---------------------------------------------------------------------!

!* Thes pointers used to arguments
  INTEGER, POINTER, DIMENSION(:,:) :: IPV       ! (NS,NV)
  INTEGER, POINTER, DIMENSION(:,:) :: IPVISO    ! (NS,NV)
  INTEGER, POINTER, DIMENSION(:)   :: ISOEFACE  ! (NV)
  INTEGER, POINTER, DIMENSION(:)   :: NIPV      ! (NS)
  INTEGER, POINTER, DIMENSION(:)   :: NIPVISO   ! (NS)
  INTEGER, POINTER                 :: NISO
  INTEGER, POINTER                 :: NTP
  INTEGER, POINTER                 :: NTS
  REAL   , POINTER, DIMENSION(:)   :: PHI       ! (NV)
  REAL   , POINTER                 :: PHIISO
  REAL   , POINTER, DIMENSION(:,:) :: VERTISO   ! (NV,3)
  REAL   , POINTER, DIMENSION(:,:) :: VERTP     ! (NV,3)
  INTEGER, POINTER, DIMENSION(:)   :: IPIA0
  INTEGER, POINTER, DIMENSION(:)   :: IPIA1

!* These are local variables
  INTEGER IA(NV),ICONTN,ICONTP,IP,IP0,IP1,IS,IV
!---------------------------------------------------------------------!

!* Transfer the pointers before all
  NIPV     => POLYHEDRON % FACES_N_NODES    ! original: NIPV
  IPV      => POLYHEDRON % FACES_N          ! original: IPV
  NTP      => POLYHEDRON % N_NODES          ! original: NTP
  NTS      => POLYHEDRON % N_FACES          ! original: NTS
  VERTP    => POLYHEDRON % NODES_XYZ        ! original: VERTP
  PHI      => POLYHEDRON % PHI
  PHIISO   => POLYHEDRON % PHI_ISO
  IPVISO   => ISO_POLYGONS % POLYS_V        ! original: IPVISO
  ISOEFACE => ISO_POLYGONS % FACE_INDEX     ! original: ISOEFACE
  NIPVISO  => ISO_POLYGONS % POLYS_N_VERTS  ! original: NIPVISO
  NISO     => ISO_POLYGONS % N_POLYS        ! original: NISO
  VERTISO  => ISO_POLYGONS % VERTS_XYZ      ! original: VERTISO
  IPIA0    => ISO_POLYGONS % B_NODE_1
  IPIA1    => ISO_POLYGONS % B_NODE_2

!* Continue as it used to be
  ICONTP=0
  ICONTN=0
  NISO=0
  DO IP=1,NTP
     IA(IP)=-1
  END DO
!* Tag the vertices of the polyhedron
  DO IS=1,NTS
     DO IV=1,NIPV(IS)
        IP=IPV(IS,IV)
        IF(IA(IP).EQ.(-1)) THEN
           IF(PHI(IP).GT.PHIISO) THEN
              IA(IP)=1
              ICONTP=ICONTP+1
           ELSE
              IA(IP)=0
              ICONTN=ICONTN+1
           END IF
        END IF
     END DO
  END DO
  IF(ICONTP.NE.0.AND.ICONTN.NE.0) THEN
!* Insert and arrange the iso-vertices
     CALL ISOAP % ISOPOL(IA,IPIA0,IPIA1,IPV,IPVISO,ISOEFACE,NIPV,  &
          NIPVISO,NISO,NTS)
!* Iso-vertices positioning by linear interpolation
     DO IS=1,NISO
        DO IV=1,NIPVISO(IS)
           IP=IPVISO(IS,IV)
           IP0=IPIA0(IP)
           IP1=IPIA1(IP)
           VERTISO(IP,1)=VERTP(IP0,1)+(PHIISO-PHI(IP0))*(       &
                VERTP(IP1,1)-VERTP(IP0,1))/(PHI(IP1)-PHI(IP0))
           VERTISO(IP,2)=VERTP(IP0,2)+(PHIISO-PHI(IP0))*(       &
                VERTP(IP1,2)-VERTP(IP0,2))/(PHI(IP1)-PHI(IP0))
           VERTISO(IP,3)=VERTP(IP0,3)+(PHIISO-PHI(IP0))*(       &
                VERTP(IP1,3)-VERTP(IP0,3))/(PHI(IP1)-PHI(IP0))
        END DO
     END DO
  END IF      
  RETURN
  END
