!=============================================================================!
  module Cgns_Mod
!-----------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
  use Sort_Mod
!-----------------------------------------------------------------------------!
  implicit none
!-----------------------------------------------------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!=============================================================================!

  !----------!
  !   File   ! -> contains base
  !----------!
    !
    !----------!
    !   Base   ! -> contains blocks
    !----------!
      !
      !------------!
      !   Blocks   ! -> contains coordinates, sections, bnd_conds, solution
      !------------!
        !
        !-----------------!
        !   Coordinates   !
        !-----------------!
        !
        !----------------------!
        !   Element sections   !
        !----------------------!
        !
        !----------------------!
        !   Solution section   !
        !----------------------!
          !
          !------------------!
          !   Field section  !
          !------------------!

  ! TO DO : 
  ! 1) Pack monitoring points in ConvergenceHistory_t
  ! 2) Salome test: mesh is converted correctly, but WallDistance is wrong

  ! File
  integer           :: file_id
  character(len=80) :: file_name
  integer           :: file_mode
  logical           :: verbose = .false.

  ! Solution section
  type Cgns_Solution_Type
    character(len=80) :: name
    integer           :: sol_type
  end type

  ! Element section
  type Cgns_Section_Type
    character(len=80)    :: name
    integer              :: cell_type
    integer              :: first_cell
    integer              :: last_cell
    integer              :: parent_flag
  end type

  ! Blocks
  type Cgns_Block_Type
    character(len=80)                       :: name
    integer                                 :: type
    integer                                 :: mesh_info(3)
    integer                                 :: n_sects
    type(Cgns_Section_Type), allocatable    :: section(:)
    integer                                 :: n_coords
    character(len=80)                       :: coord_name(3)
    integer                                 :: n_solutions
    type(Cgns_Solution_Type), allocatable   :: solution(:)
  end type

  ! Base
  integer :: n_bases
  type Cgns_Base_Type
    character(len=80)                  :: name
    integer                            :: cell_dim
    integer                            :: phys_dim
    integer                            :: n_blocks
    type(Cgns_Block_Type), allocatable :: block(:)
  end type
  type(Cgns_Base_Type), allocatable :: cgns_base(:)

  ! Some global counters (this is a bit ugly)
  integer :: cnt_nodes
  integer :: cnt_cells
  integer :: cnt_blocks     ! probably not needed
  integer :: cnt_bnd_cells
  integer :: pos_of_last_3d_cell

  ! Cells (3d)
  integer :: cnt_hex
  integer :: cnt_pyr
  integer :: cnt_wed
  integer :: cnt_tet
  integer :: cnt_mix

  ! Boundary condition cells (2d)
  integer :: cnt_bnd_qua
  integer :: cnt_bnd_tri

  ! Interface cells (2d)
  integer              :: cnt_int_qua
  integer              :: cnt_int_tri
  integer              :: cnt_int
  character(len=80)    :: interface_names(1024)
  integer, allocatable :: interface_cells(:,:,:,:)
  integer              :: cnt_int_cells

  ! Block-wise counter of boundary cells
  integer           :: cnt_block_bnd_cells  ! probably not needed
  integer           :: cnt_bnd_cond
  character(len=80) :: bnd_cond_names(1024)

  ! if an actual grid was written, further saves have just a link to that grid
  logical           :: mesh_written = .false.
  logical           :: permanent_fields_written = .false. ! like "WallDistance"
  character(len=80) :: file_with_mesh

  contains

  include 'Cgns_Mod/Initialize_Counters.f90'

  include 'Cgns_Mod/Write_Link_To_Mesh_In_File.f90'
  include 'Cgns_Mod/Write_Link_To_Field.f90'
  include 'Cgns_Mod/Write_Dimensions_Info.f90'

  ! Par only
  include 'Cgns_Mod/Parallel/Get_Arrays_Dimensions.f90'
  include 'Cgns_Mod/Parallel/Open_File.f90'
  include 'Cgns_Mod/Parallel/Close_File.f90'
  include 'Cgns_Mod/Parallel/Write_Base_Info.f90'
  include 'Cgns_Mod/Parallel/Write_Block_Info.f90'
  include 'Cgns_Mod/Parallel/Write_Coordinate_Array.f90'
  include 'Cgns_Mod/Parallel/Write_Section_Connections.f90'
  include 'Cgns_Mod/Parallel/Write_Solution_Info.f90'
  include 'Cgns_Mod/Parallel/Write_Field.f90'

  end module
