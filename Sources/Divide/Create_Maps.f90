!==============================================================================!
  subroutine Create_Maps(grid)
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod, only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, s, c1, c2, sub, n_cells_sub, n_faces_sub,  &
                          n_buf_sub, n_bnd_cells_sub, NCSsub
  character(len=80)    :: name_map
  integer, allocatable :: global_cell_ins(:)
  integer, allocatable :: global_cell_bnd(:)
!==============================================================================!

  allocate (global_cell_ins( grid % n_cells));         global_cell_ins = 0
  allocate (global_cell_bnd(-grid % n_bnd_cells:-1));  global_cell_bnd = 0

  !-------------------------------!
  !                               !
  !   Browse through subdomains   !
  !                               !
  !-------------------------------!
  do sub = 1, maxval(grid % comm % proces(:))

    call Name_File(sub, name_map, '.map')
    open(9, file=name_map)
    print *, '# Creating files: ', trim(name_map)

    ! Cells
    n_cells_sub = 0     ! number of cells in subdomain
    do c = 1, grid % n_cells
      if(grid % comm % proces(c) .eq. sub) then
        n_cells_sub = n_cells_sub + 1     ! increase the number of cells in sub.
        global_cell_ins(n_cells_sub) = c  ! map to global cell number
      end if
    end do

    ! Faces & real boundary cells
    n_faces_sub     = 0  ! number of sides in subdomain
    n_buf_sub       = 0  ! number of buffer faces (and cells)
    n_bnd_cells_sub = 0  ! number of real boundary cells in subdomain
    NCSsub = 0

    ! Faces step 2/3: on the boundaries + bundary cells
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)  
      c2 = grid % faces_c(2,s) 
      if(c2 < 0) then
        if( grid % comm % proces(c1) .eq. sub )  then
          n_bnd_cells_sub = n_bnd_cells_sub + 1  ! increase n. of bnd. cells
          global_cell_bnd(-n_bnd_cells_sub) = c2  ! map to global cell number
        end if
      end if 
    end do

    !-----------------------!
    !   Save cell mapping   !
    !-----------------------!

    ! First line are the number of boundary cells and cells
    write(9, '(4i9)') n_cells_sub, n_bnd_cells_sub

    ! Extents are followed by mapping of the cells inside ...
    do c = 1, n_cells_sub
      write(9, '(i9)') global_cell_ins(c)
    end do

    ! ... followed by the cells on the boundary ...
    do c = -n_bnd_cells_sub, -1
      write(9, '(i9)') global_cell_bnd(c)
    end do

  end do   ! through subdomains

  close(9)

  end subroutine
