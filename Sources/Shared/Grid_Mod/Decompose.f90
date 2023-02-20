!==============================================================================!
  subroutine Decompose(Grid, n_parts)
!------------------------------------------------------------------------------!
!   Coarsens the grid with METIS library.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: n_parts
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, s
  integer, allocatable :: part(:)          ! result of partitioning
!==============================================================================!

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  allocate(part(Grid % n_cells))

  !----------------------------------------------------------------------------!
  !   Prepare METIS and Exectute a call to METIS graph partitioning function   !
  !----------------------------------------------------------------------------!
  call Metis % Create_Metis(Grid % faces_c, n_parts)

  Metis % options = -1                       ! Initialize all to default
  Metis % options(METIS_OPTION_DBGLVL) = 0

  call Metis_PartGraphRecursive(     &  ! rank  intent  type     METIS
              Metis % n_verts,       &  !   1.  (in),   int      nvtxs
              Metis % n_constrains,  &  !   2.  (in),   int      ncon
              Metis % row,           &  !   3.  (in),   int(:)   xadj
              Metis % col,           &  !   4.  (in),   int(:)   adjncy
              Metis % vert_weights,  &  !   5.  (in),   int(:)   vwgt
              Metis % vert_data,     &  !   6.  (in),   int(:)   vsize
              Metis % edge_weights,  &  !   7.  (in),   int(:)   adjwgt
              n_parts,               &  !   8.  (in),   int(:)   nparts
              Metis % part_weight,   &  !   9.  (in),   real(:)  tpwgts
              Metis % imbalance,     &  !  10.  (in),   real(:)  ubvec
              Metis % options,       &  !  11.  (in),   int(:)   options
              Metis % return_val,    &  !  12.  (out)   int(:)   objval
              part(:))                  !  13.  (out)   int(:)   part

  part(:) = part(:) + 1  ! +1, METIS works from zero

  !-----------------------------------------------------!
  !   Save the result from the call to METIS function   !
  !-----------------------------------------------------!
  do c = 1, Grid % n_cells
    Grid % Comm % cell_proc(c) = part(c)
  end do

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) then
      Grid % Comm % cell_proc(c2) = Grid % Comm % cell_proc(c1)
    end if
  end do

  end subroutine
