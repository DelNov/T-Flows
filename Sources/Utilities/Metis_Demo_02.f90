!==============================================================================!
  program Test_Metis_02
!------------------------------------------------------------------------------!
!   Second demonstration program which calls Metis library.                    !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Iso_C_Binding
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
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
!-----------------------------------[Locals]-----------------------------------!
  integer     :: nvtxs,              &  ! number of vertices
                 ncons,              &  ! number of constrains 
                 nparts = 2             ! requested number of partitions
  integer        mopts(41)
  type(c_ptr) :: vwgt  =c_null_ptr,  &  ! weights of the vertices
                 vsize =c_null_ptr,  &  ! size of the vertices
                 adjwgt=c_null_ptr,  &  ! weights of the edges
                 objval=c_null_ptr      ! stores edgecut or comm. volume
  real,    allocatable :: tpwgts(:)          ! desired weight for each partition
  real,    allocatable :: ubvec (:)
  integer, allocatable :: xadj  (:),      &  ! variabes for ...
                          adjncy(:),      &  ! ... compressed row storage
                          part  (:)          ! partition of the grid
!==============================================================================!

  nvtxs = 15
  allocate(xadj  (nvtxs+1))
  allocate(adjncy(22 * 2))
  allocate(part  (nvtxs))

  ncons = 1
  allocate(tpwgts(nparts*ncons))
  allocate(ubvec (ncons))

  xadj = (/ 0,  2,  5,  8, 11, 13, 16, 20,  &
           24, 28, 31, 33, 36, 39, 42, 44/)

  adjncy=(/ 1,  5,  0,  2,  6,  1,  3,  7,  2,  4,  8,  &
            3,  9,  0,  6, 10,  1,  5,  7, 11,  2,  6,  &
            8, 12,  3,  7,  9, 13,  4,  8, 14,  5, 11,  &
            6, 10, 12,  7, 11, 13,  8, 12, 14,  9, 13/)
  tpwgts(:) = 1.0 / real(nparts)
  ubvec(:)  = 1.01

  call METIS_SetDefaultOptions(mopts)
  mopts(METIS_OPTION_DBGLVL)    = 1

  ! Call METIS function
  call METIS_PartGraphRecursive(nvtxs,     &  ! (in), int
                                ncons,     &  ! (in), int
                                xadj,      &  ! (in), int(:)
                                adjncy,    &  ! (in), int(:)
                                vwgt,      &  ! (in), int(:)
                                vsize,     &  ! (in), int(:)
                                adjwgt,    &  ! (in), int(:)
                                nparts,    &  ! (in), int(:)
                                tpwgts,    &  ! (in), real(:)
                                ubvec,     &  ! (in), real(:)
                                mopts,     &  ! (in), int(:)
                                objval,    &  ! (out) int(:)
                                part)         ! (out) int(:)
 
  print *, part

  end program
