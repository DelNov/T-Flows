!==============================================================================!
  subroutine Call_Metis(Metis, n_parts, part)
!------------------------------------------------------------------------------!
!>  Facilitates the partitioning of grid into a specified number of parts using
!>  the METIS library. This subroutine abstracts the complexity of calling
!>  METIS functions by providing a simplified interface
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Metis_Type),     intent(inout) :: Metis    !! parent class
  integer,               intent(in)    :: n_parts  !! number of partitions
  integer, dimension(:), intent(out)   :: part     !! result of partitioning
!==============================================================================!

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

  end subroutine
