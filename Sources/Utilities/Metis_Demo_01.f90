!==============================================================================!
  program Test_Metis_01
!------------------------------------------------------------------------------!
!   First demonstration program which calls Metis library.                     !
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Locals]------------------------------------!
  integer              :: nels   = 2,  &  ! number of elements
                          nnds   = 6,  &  ! number of nodes
                          npel   = 4,  &  ! nodes per element
                          nparts = 2,  &  ! number of desired partitions
                          objval          ! returns something after METIS
  integer, allocatable :: eptr(:),   &
                          nodes(:),  &
                          epart(:),  &
                          npart(:)              
  integer, pointer     :: vwgt   => null(),  &
                          vsize  => null()
  real, pointer        :: tpwgts => null()
  integer              :: mopts(41)       ! options passed to METIS 

  integer, parameter :: METIS_OPTION_PTYPE     =  1
  integer, parameter :: METIS_OPTION_OBJTYPE   =  2
  integer, parameter :: METIS_OPTION_CTYPE     =  3
  integer, parameter :: METIS_OPTION_IPTYPE    =  4
  integer, parameter :: METIS_OPTION_RTYPE     =  5
  integer, parameter :: METIS_OPTION_DBGLVL    =  6
  integer, parameter :: METIS_OPTION_NITER     =  7
  integer, parameter :: METIS_OPTION_NCUTS     =  8
  integer, parameter :: METIS_OPTION_SEED      =  9
  integer, parameter :: METIS_OPTION_NO2HOP    = 10
  integer, parameter :: METIS_OPTION_MINCONN   = 11
  integer, parameter :: METIS_OPTION_CONTIG    = 12
  integer, parameter :: METIS_OPTION_COMPRESS  = 13
  integer, parameter :: METIS_OPTION_CCORDER   = 14
  integer, parameter :: METIS_OPTION_PFACTOR   = 15
  integer, parameter :: METIS_OPTION_NSEPS     = 16
  integer, parameter :: METIS_OPTION_UFACTOR   = 17
  integer, parameter :: METIS_OPTION_NUMBERING = 18
  integer, parameter :: METIS_OPTION_HELP      = 19
  integer, parameter :: METIS_OPTION_TPWGTS    = 20
  integer, parameter :: METIS_OPTION_NCOMMON   = 21
  integer, parameter :: METIS_OPTION_NOOUTPUT  = 22
  integer, parameter :: METIS_OPTION_BALANCE   = 23
  integer, parameter :: METIS_OPTION_GTYPE     = 24
  integer, parameter :: METIS_OPTION_UBVEC     = 25
!------------------------------------------------------------------------------!

  allocate(eptr(nels+1))
  allocate(nodes(nels*npel))
  allocate(epart(nels))
  allocate(npart(nnds))

  ! Initialize data to be sent to Metis
  eptr  = (/0, 4, 8/)
  nodes = (/1, 2, 3, 4, 2, 5, 6, 3/)  ! element 1 has nodes 1 2 3 4
                                      ! element 2 has nodes 2 5 6 3

  ! Set options to default values
  mopts = -1
  mopts(METIS_OPTION_NUMBERING) = 1
  mopts(METIS_OPTION_DBGLVL)    = 1

  ! Call Metis
  call METIS_PartMeshNodal(nels,    &
                           nnds,    &
                           eptr,    &
                           nodes,   &
                           vwgt,    &
                           vsize,   &
                           nparts,  &
                           tpwgts,  &
                           mopts,   &
                           objval,  &
                           epart,   &
                           npart) 

  print *, '# Node partitioning: ', npart  
  print *, '# Elem partitioning: ', epart

  end program
