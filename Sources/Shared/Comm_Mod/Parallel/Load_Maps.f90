!==============================================================================!
  subroutine Comm_Mod_Load_Maps(grid)
!------------------------------------------------------------------------------!
!   Reads: name.map file                                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c
  character(len=80) :: name_in
!==============================================================================!
!   There is an issue with this procedure, but it's more related to MPI/IO     !
!   functions than T-Flows.  In cases a subdomain has no physical boundary     !
!   cells, variable "nb_s" turns out to be zero.  This, per se, should not     !
!   be an issue if MPI/IO functions could handle calls to:                     !
!   "Mpi_Type_Create_Indexed_Block(...)" and "Mpi_File_Write(...)" with zero   !
!   length.  But they don't.  Therefore, I avoid allocation with zero size     !
!   (max(nb_s,1)) here and creation of new types with zero size in             !
!   "Comm_Mod_Create_New_Types".  It is a bit of a dirty trick :-(             !
!------------------------------------------------------------------------------!

  !------------------------------------------------------------------------!
  !                                                                        !
  !   For run with one processor, no needd to read the map, just form it   !
  !                                                                        !
  !------------------------------------------------------------------------!
  if(n_proc < 2) then

    nc_s  = grid % n_cells
    nb_s  = grid % n_bnd_cells
    nc_t  = nc_s
    nb_t  = nb_s

    !-------------------------------------!
    !   Global cell numbers for T-Flows   !
    !-------------------------------------!
    allocate(grid % comm % cell_glo(-nb_s:nc_s))

    ! Fill up global cell numbers
    do c = -nb_t, nc_t
      grid % comm % cell_glo(c) = c
    end do

    !-----------------------------------------!
    !   Global cell numbers for MPI mapping   !
    !-----------------------------------------!
    allocate(grid % comm % cell_map    (nc_s))  
    allocate(grid % comm % bnd_cell_map(nb_s))

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, nc_t
      grid % comm % cell_map(c) = c - 1
    end do
  
    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, nb_t
      grid % comm % bnd_cell_map(c) = c - 1
    end do

  !-------------------------------------------------!
  !                                                 !
  !   For parallel runs, you need to read the map   !
  !                                                 !
  !-------------------------------------------------!
  else

    call File_Mod_Set_Name(name_in, processor=this_proc, extension='.map')
    open(9, file=name_in)
    if(this_proc < 2) print *, '# Now reading the file:', name_in

    ! Read map sizes
    read(9, '(4i9)') nc_s, nb_s

    nc_t  = nc_s
    nb_t  = nb_s
    call Comm_Mod_Global_Sum_Int(nc_t)
    call Comm_Mod_Global_Sum_Int(nb_t)

    !-------------------------------------!
    !   Global cell numbers for T-Flows   !
    !-------------------------------------!
    allocate(grid % comm % cell_glo(-grid % n_bnd_cells : grid % n_cells))
    grid % comm % cell_glo(:) = 0

    !-----------------------------------------!
    !   Global cell numbers for MPI mapping   !
    !-----------------------------------------!
    allocate(grid % comm % cell_map    (nc_s))   
    allocate(grid % comm % bnd_cell_map(max(nb_s,1)))  ! avoid zero allocation   
    grid % comm % cell_map(:)     = 0       
    grid % comm % bnd_cell_map(:) = 0

    !-------------------!
    !   Read cell map   !
    !-------------------!
    do c = 1, nc_s
      read(9, '(i9)') grid % comm % cell_glo(c)

      ! Take cell mapping to be the same as global cell numbers but start from 0
      grid % comm % cell_map(c) = grid % comm % cell_glo(c) - 1
    end do

    !----------------------------!
    !   Read boundary cell map   !
    !----------------------------!
    do c = -nb_s, -1
      read(9, '(i9)') grid % comm % cell_glo(c)

      ! Correct boundary cell mapping.  
      ! - First it is in positive range, so insted of -nb_s to -1, it goes from 
      !   1 to nb_s.  (Therefore the "c+nb_s+1")
      ! - Second, mapping must be positive and start from zero.  (The "+ nb_t")
      grid % comm % bnd_cell_map(c+nb_s+1) = grid % comm % cell_glo(c) + nb_t
    end do
    
    !---------------------------------------------!
    !   Refresh buffers for global cell numbers   !
    !---------------------------------------------!
    call Comm_Mod_Exchange_Int(grid, grid % comm % cell_glo)

    close(9)

  end if

  end subroutine
