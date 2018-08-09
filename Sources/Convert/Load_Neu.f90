!==============================================================================!
  subroutine Load_Neu(grid) 
!------------------------------------------------------------------------------!
!   Reads the Fluents (Gambits) neutral file format.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use gen_mod 
  use Grid_Mod
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  character(len=130)  :: name_in
  integer             :: i, j, n_blocks, n_bnd_sect, dum1, dum2
  integer,allocatable :: temp(:)
  integer             :: c, dir
!==============================================================================!

  name_in = problem_name
  name_in(len_trim(problem_name)+1:len_trim(problem_name)+4) = '.neu'

  open(9,file=name_in)
  print *, '# Reading the file: ', trim(name_in)

  ! Skip first 6 lines
  do i = 1, 6
    call Tokenizer_Mod_Read_Line(9)
  end do 

  ! Read the line which contains usefull information  
  call Tokenizer_Mod_Read_Line(9)

  read(line % tokens(1),*) grid % n_nodes  
  read(line % tokens(2),*) grid % n_cells
  read(line % tokens(3),*) n_blocks
  read(line % tokens(4),*) n_bnd_sect

  print '(a38,i9)', '# Total number of nodes:             ', grid % n_nodes
  print '(a38,i9)', '# Total number of cells:             ', grid % n_cells
  print '(a38,i9)', '# Total number of blocks:            ', n_blocks
  print '(a38,i9)', '# Total number of boundary sections: ', n_bnd_sect

  !------------------------------!
  !   Count the boundary cells   !
  !------------------------------!
  grid % n_bnd_cells = 0
  do 
    call Tokenizer_Mod_Read_Line(9)
    if( line % tokens(1) .eq. 'BOUNDARY' ) then
      do j = 1, n_bnd_sect
        if(j>1) call Tokenizer_Mod_Read_Line(9) ! BOUNDARY CONDITIONS
        call Tokenizer_Mod_Read_Line(9)
        read(line % tokens(3),*) dum1  
        grid % n_bnd_cells = grid % n_bnd_cells + dum1 
        do i = 1, dum1
          read(9,*) c, dum2, dir
        end do
        call Tokenizer_Mod_Read_Line(9)         ! ENDOFSECTION
      end do
      print '(a38,i9)', '# Total number of boundary cells:    ',  &
            grid % n_bnd_cells
      go to 1
    end if
  end do 

1 rewind(9)

  !------------------------!
  !   Skip first 7 lines   !
  !------------------------!
  do i = 1, 7
    call Tokenizer_Mod_Read_Line(9)
  end do 

  !--------------------------------------------!
  !                                            !
  !   Allocate memory for Grid_Mod variables   !
  !                                            !
  !--------------------------------------------!
  call Allocate_Memory(grid)

  allocate(temp(grid % n_cells)); temp=0

  ! Skip one line 
  call Tokenizer_Mod_Read_Line(9)

  !--------------------------------!
  !   Read the nodal coordinates   !
  !--------------------------------!
  call Tokenizer_Mod_Read_Line(9)          ! NODAL COORDINATES
  do i = 1, grid % n_nodes
    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(2),*) grid % xn(i) 
    read(line % tokens(3),*) grid % yn(i)
    read(line % tokens(4),*) grid % zn(i)
  end do
  call Tokenizer_Mod_Read_Line(9)          ! ENDOFSECTION

  !-----------------------------!
  !   Read nodes of each cell   !
  !-----------------------------!
  call Tokenizer_Mod_Read_Line(9)          ! ELEMENTS/CELLS
  do i = 1, grid % n_cells
    read(9,'(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))') dum1, dum2, &
           grid % cells_n_nodes(i),                           &
          (grid % cells_n(j,i), j = 1,grid % cells_n_nodes(i))
  end do
  call Tokenizer_Mod_Read_Line(9)          ! ENDOFSECTION

  !---------------------------------!
  !   Read block (material?) data   !
  !---------------------------------!
  allocate(grid % materials(n_blocks))
  if(n_blocks .gt. 1) then
    print *, '# Multiple materials from .neu file not been implemented yet!'
    print *, '# Exiting!'
  end if

  grid % n_materials = 1
  grid % materials(1) % name = "AIR"

  do j = 1, n_blocks
    call Tokenizer_Mod_Read_Line(9)        ! ELEMENT GROUP
    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(4),'(I10)') dum1  
    call Tokenizer_Mod_Read_Line(9)        ! block*
    call Tokenizer_Mod_Read_Line(9)        ! 0
    read(9,'(10I8)') (temp(i), i = 1, dum1)
    do i = 1, dum1
      grid % material(temp(i)) = j
    end do
    call Tokenizer_Mod_Read_Line(9)        ! ENDOFSECTION
  end do

  !-------------------------!
  !   Boundary conditions   !
  !-------------------------!
  grid % n_bnd_cond = n_bnd_sect
  allocate(grid % bnd_cond % name(n_bnd_sect))

  do j = 1, n_bnd_sect
    call Tokenizer_Mod_Read_Line(9)        ! BOUNDARY CONDITIONS
    call Tokenizer_Mod_Read_Line(9)
    call To_Upper_Case(  line % tokens(1)  )
    grid % bnd_cond % name(j) = line % tokens(1)
    read(line % tokens(3),'(I8)') dum1  
    do i = 1, dum1
      read(9,*) c, dum2, dir
      grid % cells_bnd_color(dir,c) = j
    end do
    call Tokenizer_Mod_Read_Line(9)        ! ENDOFSECTION
  end do

  !------------------------------------!
  !   Pring boundary conditions info   !
  !------------------------------------!
  call Grid_Mod_Print_Bnd_Cond_List(grid)

  close(9)

  end subroutine
