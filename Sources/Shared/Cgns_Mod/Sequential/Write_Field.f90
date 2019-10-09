!==============================================================================!
  subroutine Cgns_Mod_Write_Field(base, block, solution, field, grid, &
    input_array, input_name)
!------------------------------------------------------------------------------!
!   Writes field to solution node and sets its field_id  [sequential vesion]   !
!------------------------------------------------------------------------------!
!   Array structures in current function are strictly followings:              !
!                                                                              !
!   Cell type:    |      HEXA_8      |     PENTA_6      |       PYRA_5     |...!
!                                                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base, block, solution, field
  type(Grid_Type)      :: grid
  real                 :: input_array(grid % n_cells)
  character(len=*)     :: input_name
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id     ! base index number
  integer              :: block_id    ! block index number
  integer              :: solution_id ! solution index
  integer              :: field_id    ! field index
  character(len=80)    :: field_name  ! name of the FlowSolution_t node
  integer              :: sect_id
  integer              :: cnt            ! cells of sect_id
  real                 :: field_array(grid % n_cells) ! field array
  integer              :: i, j, k, c
  integer              :: error
!==============================================================================!

  ! Set input parameters
  base_id       = base
  block_id      = block
  solution_id   = solution
  field_id      = field

  field_name = trim(input_name)

  !------------------------------------------------!
  !   Mapping 1:nc -> Connection structure above   !
  !------------------------------------------------!

  ! Find first and last cells of sect_id
  do sect_id = 1, 4
    ! cells of sect_id
    if (sect_id.eq.1) cnt = cnt_hex
    if (sect_id.eq.2) cnt = cnt_wed
    if (sect_id.eq.3) cnt = cnt_pyr
    if (sect_id.eq.4) cnt = cnt_tet

    if (cnt.ne.0) then

      i = 1 + cnt_cells
      j = i + cnt - 1

      ! copy input array to field_array
      k = 1
      do c = 1, grid % n_cells - grid % comm % n_buff_cells
        if     (sect_id.eq.1 .and. grid % cells_n_nodes(c).eq.8) then
          field_array(k) = input_array(c)
          k = k + 1
        elseif (sect_id.eq.2 .and. grid % cells_n_nodes(c).eq.6) then
          field_array(k) = input_array(c)
          k = k + 1
        elseif (sect_id.eq.3 .and. grid % cells_n_nodes(c).eq.5) then
          field_array(k) = input_array(c)
          k = k + 1
        elseif (sect_id.eq.4 .and. grid % cells_n_nodes(c).eq.4) then
          field_array(k) = input_array(c)
          k = k + 1
        end if
      end do

      !----------------------------------------------!
      !   Writing cells assocciated with sect_id   !
      !----------------------------------------------!

      ! Add field to FlowSolution_t node for sect_id
      call Cg_Field_Partial_Write_F(file_id,      & !(in )
                                    base_id,      & !(in )
                                    block_id,     & !(in )
                                    solution_id,  & !(in )
                                    RealDouble,   & !(in )
                                    field_name,   & !(in )
                                    i,            & !(in )
                                    j,            & !(in )
                                    field_array,  & !(in )
                                    field_id,     & !(out)
                                    error)          !(out)

      if (error .ne. 0) then
        print *, '# Failed to write field', trim(field_name)
        call Cg_Error_Exit_F()
      endif

      cnt_cells = cnt_cells + cnt

      ! Print some info
      if(verbose ) then
        print *, '#           ---------------------------------'
        print *, '#           Field name: ',   field_name
        print *, '#           Field idx:    ', field_id
        print *, '#           ---------------------------------'
        print *, '#           First cell:', i
        print *, '#           Last cell: ', j
      end if

    end if

  end do

  cnt_cells = 0

  end subroutine
