!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for MPI functionality.                                              !
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
  include 'mpif.h'
!==============================================================================!

  ! Parameters for size of typical variables in bytes
  integer, parameter :: SIZE_INT  = 8
  integer, parameter :: SIZE_REAL = 8

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  ! These names are ugly but mean number of buffer boundaries start and end
  integer, allocatable :: nbb_s(:), nbb_e(:)
  integer, allocatable :: buffer_index(:)

  ! Variables which follow are for backup saving to single file
  integer                      :: nc_s   ! number of cells in subdomain
  integer                      :: nb_s   ! number of bundary cells in subdoima
  integer                      :: nf_s   ! number of faces in subdomain
  integer                      :: nbf_s  ! number of buffer faces in subdomain
  integer                      :: nc_t   ! total number of cells 
  integer                      :: nb_t   ! total number of bundary cells
  integer                      :: nf_t   ! total number of faces
  integer(kind=4), allocatable :: cell_map(:)
  integer(kind=4), allocatable :: bnd_cell_map(:)
  integer(kind=4), allocatable :: face_map(:)
  integer(kind=4), allocatable :: face_ord(:)
  real, allocatable            :: face_val(:)      ! face values
  integer(kind=4), allocatable :: buf_face_map(:)  ! buffer face map
  integer(kind=4), allocatable :: buf_cell_map(:)  ! buffer cell map
  integer(kind=4), allocatable :: buf_face_ord(:)  ! buffer face order
  real, allocatable            :: buf_face_val(:)  ! buffer face values
  real, allocatable            :: buf_face_sgn(:)  ! buffer face sign  
  integer                      :: cell_map_type
  integer                      :: bnd_cell_map_type
  integer                      :: face_map_type
  integer                      :: buf_face_map_type

  contains

  include 'Comm_Mod/Parallel/Allocate_Memory.f90'
  include 'Comm_Mod/Parallel/Close_File.f90'
  include 'Comm_Mod/Parallel/Create_New_Types.f90'
  include 'Comm_Mod/Parallel/End.f90'
  include 'Comm_Mod/Parallel/Exchange.f90'
  include 'Comm_Mod/Parallel/Global_Max_Real.f90'
  include 'Comm_Mod/Parallel/Global_Min_Real.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Real.f90'
  include 'Comm_Mod/Parallel/Load_Buffers.f90'
  include 'Comm_Mod/Parallel/Load_Maps.f90'
  include 'Comm_Mod/Parallel/Open_File_Read.f90'
  include 'Comm_Mod/Parallel/Open_File_Write.f90'
  include 'Comm_Mod/Parallel/Read_Int.f90'
  include 'Comm_Mod/Parallel/Read_Bnd_Real.f90'
  include 'Comm_Mod/Parallel/Read_Cell_Real.f90'
  include 'Comm_Mod/Parallel/Read_Face_Real.f90'
  include 'Comm_Mod/Parallel/Read_Real.f90'
  include 'Comm_Mod/Parallel/Read_Text.f90'
  include 'Comm_Mod/Parallel/Start.f90'
  include 'Comm_Mod/Parallel/Wait.f90'
  include 'Comm_Mod/Parallel/Write_Int.f90'
  include 'Comm_Mod/Parallel/Write_Bnd_Real.f90'
  include 'Comm_Mod/Parallel/Write_Cell_Real.f90'
  include 'Comm_Mod/Parallel/Write_Face_Real.f90'
  include 'Comm_Mod/Parallel/Write_Real.f90'
  include 'Comm_Mod/Parallel/Write_Text.f90'

  end module
