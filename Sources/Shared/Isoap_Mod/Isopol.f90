!---------------------------------------------------------------------!
!                              ISOPOL                                 !
!---------------------------------------------------------------------!
!>  Original Isoap's subroutine to insert and arrange the iso-vertices,
!>  called from main Isoap function.
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
! IA       = tag value of each vertex of the polyhedron (0 or 1)      !
! IPV      = array containing the indices of the vertices of each     !
!            face of the polyhedron                                   !
! NIPV     = number of vertices of each face of the polyhedron        !
! NTS      = number of faces of the polyhedron                        !
! On return:                                                          !
!===========                                                          !
! IPIA0    = vertex index of the polihedron with IA=0 and which is in !
!            the edge containing the iso-vertex                       !
! IPIA1    = vertex index of the polihedron with IA=1 and which is in !
!            the edge containing the iso-vertex                       !
! IPVISO   = array conatining the indices of the iso-vertices of each !
!            iso-polygon                                              !
! ISOEFACE = face index of the polyhedron over which is constructed   !
!            each iso-edge                                            !
! NIPVISO  = number of iso-vertices of each iso-polygon               !
! NISO     = number of iso-polygons                                   !
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
  SUBROUTINE ISOPOL(ISOAP,     &
                    IA,        &
                    IPIA0,     &
                    IPIA1,     &
                    IPV,       &
                    IPVISO,    &
                    ISOEFACE,  &
                    NIPV,      &
                    NIPVISO,   &
                    NISO,      &
                    NTS)
!---------------------------------------------------------------------!
  IMPLICIT NONE
!---------------------------------------------------------------------!
  CLASS(ISOAP_TYPE)    :: ISOAP
  INTEGER, INTENT(IN)  :: IA(NV)
  INTEGER, INTENT(OUT) :: IPIA0(NV)
  INTEGER, INTENT(OUT) :: IPIA1(NV)
  INTEGER, INTENT(IN)  :: IPV(NS,NV)
  INTEGER, INTENT(OUT) :: IPVISO(NS,NV)
  INTEGER, INTENT(OUT) :: ISOEFACE(NV)
  INTEGER, INTENT(IN)  :: NIPV(NS)
  INTEGER, INTENT(OUT) :: NIPVISO(NS)
  INTEGER, INTENT(OUT) :: NISO
  INTEGER, INTENT(IN)  :: NTS
!---------------------------------------------------------------------!
!* Iso-vertices
  INTEGER IPVINT(NS,NV),NIPVINT(NS)
!* Working parameters
  INTEGER IE,IP,IP0I,IP0N,IP1,IP1N,IP1I,  &
       IPINI,IPISE(NV,2),IPMARK(NV),IPNEW,IS,IS1,ISCUT(NS),          &
       ISE(NS,NV),ISNEW,ITYPE,IV,IV1,IVISE(NS,NV),IVNEW,IVNEWT,      &
       NEDGE(NS),NINT,NIPNEW,NISCUT,NIV,NIVNEW
!* Avoid unused warning
  ASSOCIATE(ISOAP => ISOAP); END ASSOCIATE
!* Determination of the faces intersected by the isosurface      
  NISCUT=0
!* NEDGE(IS) = Number of intersected edges of the face IS
  DO IS=1,NTS
     NEDGE(IS)=0
     IF(NIPV(IS).GT.0) THEN
        ISCUT(IS)=0
        DO IV=1,NIPV(IS)
           IP=IPV(IS,IV)
           IV1=IV+1
           IF(IV.EQ.NIPV(IS)) IV1=1
           IP1=IPV(IS,IV1)
           IF(IA(IP).NE.IA(IP1)) THEN
              ISCUT(IS)=1
              NISCUT=NISCUT+1
              NEDGE(IS)=NEDGE(IS)+1
           END IF
        END DO
     END IF
  END DO
!* Disjoint regions may produce NISCUT=0 and both ICONTP and ICONTN \NEQ 0
  IF(NISCUT.EQ.0) THEN
     NISO=0
     RETURN
  END IF      
!* Iso-vertices insertion
  NIPNEW=0
  DO IS=1,NTS
     IF(ISCUT(IS).EQ.1) THEN
        NIV=0
        NINT=0
        DO IV=1,NIPV(IS)
           IP=IPV(IS,IV)
           IV1=IV+1
           IF(IV1.GT.NIPV(IS))IV1=1
           IP1=IPV(IS,IV1)               
           IF(IA(IP).NE.IA(IP1)) THEN
              NINT=NINT+1
              NIV=NIV+1
              IF(IA(IP).EQ.1) THEN
                 IP1I=IP
                 IP0I=IP1
                 ITYPE=2
              ELSE
                 IP1I=IP1
                 IP0I=IP
                 ITYPE=1
              END IF
              DO IS1=1,IS-1
                 DO IE=1,NEDGE(IS1)
                    IPNEW=ISE(IS1,IE)
                    IP0N=IPIA0(IPNEW)
                    IP1N=IPIA1(IPNEW)
                    IF(IP0N.EQ.IP0I.AND.IP1N.EQ.IP1I) THEN
                       ISE(IS,NINT)=IPNEW
                       IPVINT(IS,NIV)=IPNEW                           
                       IVISE(IS,IPNEW)=NIV
                       IPISE(IPNEW,ITYPE)=IS
                       GOTO 10
                    END IF
                 END DO     
              END DO        
              NIPNEW=NIPNEW+1
              IPIA0(NIPNEW)=IP0I
              IPIA1(NIPNEW)=IP1I
              IPVINT(IS,NIV)=NIPNEW
              ISE(IS,NINT)=NIPNEW
              IVISE(IS,NIPNEW)=NIV
              IPISE(NIPNEW,ITYPE)=IS
           END IF           
10            CONTINUE
        END DO
        NIPVINT(IS)=NIV
     END IF                 
  END DO                    
!* Iso-vertices arrangement
  NIVNEW=NIPNEW
  ISNEW=0
  DO IP=1,NIPNEW
     IPMARK(IP)=0
  END DO
  IVNEWT=0
  IPNEW=1
!* First point
40   CONTINUE
  IVNEW=1
  IVNEWT=IVNEWT+1
  ISNEW=ISNEW+1
  IPINI=IPNEW
  IPVISO(ISNEW,IVNEW)=IPNEW
  ISOEFACE(IPNEW)=ipise(ipnew,1)
  IPMARK(IPNEW)=1
20   CONTINUE
  IS=IPISE(IPNEW,1)
  IV=IVISE(IS,IPNEW)
  IV1=IV-1
  IF(IV1.EQ.0) IV1=NIPVINT(IS)
  IPNEW=IPVINT(IS,IV1)
  IF(IPNEW.NE.IPINI) THEN
     IVNEW=IVNEW+1
     IVNEWT=IVNEWT+1
     IPVISO(ISNEW,IVNEW)=IPNEW
     ISOEFACE(IPNEW)=ipise(ipnew,1)
     IPMARK(IPNEW)=1
     IF(IVNEWT.EQ.NIVNEW) GOTO 30
     GOTO 20
  END IF
  NIPVISO(ISNEW)=IVNEW
  DO IPNEW=2,NIPNEW
     IF(IPMARK(IPNEW).EQ.0) GOTO 40
  END DO
30   CONTINUE
  NIPVISO(ISNEW)=IVNEW
  NISO=ISNEW

  RETURN
  END
