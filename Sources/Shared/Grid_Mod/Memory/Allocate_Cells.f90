!==============================================================================!
  subroutine Allocate_Cells(Grid, nc, nb, margin)
!------------------------------------------------------------------------------!
!>  Allocates memory for cell-based data (arrays and matrices), for geometrical
!>  (xc, yc, zc, vol ...) and connectivity data (cells_n, cells_f, ...).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)              :: Grid    !! grid being expanded
  integer, intent(in)           :: nc      !! number of cells inside
  integer, intent(in)           :: nb      !! number of cells on the bounday
  integer, intent(in), optional :: margin  !! margin for allocation
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, nc_m, nb_m
!==============================================================================!

  ! Generator has growing number of cells, don't set them here
  if(PROGRAM_NAME .ne. 'Generate') then
    Grid % n_cells     = nc
    Grid % n_bnd_cells = nb
  end if

  ! If cell-based arrays are allocated and they
  ! are bigger than requested, get out of here
  if(allocated(Grid % xc)) then
    if(-nb .ge. lbound(Grid % xc, 1) .and. nc .le. ubound(Grid % xc, 1)) then
      return
    end if
  end if

  ! Process the margin if specified
  nc_m = nc
  nb_m = nb
  if(present(margin)) then
    Assert(margin .ge. 0)
    nc_m = nc + margin
    nb_m = nb + margin
  end if

  ! I can't figure out how to print something meaningful in parallel here
  if(Sequential_Run()) then
    print '(a,2i9)', ' # Expanding memory for cells to size: ', nc, nb
  end if

  ! Allocate cell center coordinates and initialize to zero
  call Enlarge % Array_Real(Grid % xc, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % yc, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % zc, i=(/-nb_m,nc_m/))

  ! Memory for cells' volumes, delta and wall distance
  call Enlarge % Array_Real(Grid % vol,            i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % wall_dist,      i=(/-nb_m,nc_m/))
  call Enlarge % Array_Log (Grid % cell_near_wall, i=(/-nb_m,nc_m/))

  ! Memory for cells' inertia tensors
  call Enlarge % Array_Real(Grid % ixx, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % iyy, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % izz, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % ixy, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % ixz, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Real(Grid % iyz, i=(/-nb_m,nc_m/))

  ! Allocate as litle as possible
  call Enlarge % Matrix_Int(Grid % cells_n, i=(/1,4/), j=(/-nb_m,nc_m/))
  call Enlarge % Matrix_Int(Grid % cells_f, i=(/1,4/), j=(/-nb_m,nc_m/))
  if(PROGRAM_NAME      .eq. 'Generate' .or.  &
     PROGRAM_NAME      .eq. 'Convert'  .or.  &
     PROGRAM_NAME(1:7) .eq. 'Process') then
    call Enlarge % Matrix_Int(Grid % cells_c, i=(/1,4/), j=(/-nb_m,nc_m/))
  end if

  ! This is used only in Swarm_Mod for bouncing particles
  if(PROGRAM_NAME(1:7) .eq. 'Process') then
    if(nb_m .gt. 0) then  ! this check is not crazy, a (sub)domain might
                          ! have no boundary cells, all could be periodic
                          ! or simply not touching any boundary regions
      call Enlarge % Array_Int(Grid % cells_bnd_face, i=(/-nb_m,-1/))
    end if
  end if

  ! Number of nodes, faces and cells at each cell
  ! (Actually, cells_n_faces and cells_n_cells should be the same)
  call Enlarge % Array_Int(Grid % cells_n_nodes, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Int(Grid % cells_n_faces, i=(/-nb_m,nc_m/))
  if(PROGRAM_NAME      .eq. 'Generate' .or.  &
     PROGRAM_NAME      .eq. 'Convert'  .or.  &
     PROGRAM_NAME(1:7) .eq. 'Process') then
    call Enlarge % Array_Int(Grid % cells_n_cells, i=(/-nb_m,nc_m/))
  end if

  ! Boundary condition region in a given direction
  ! (These go up to 6 because they are needed for
  !  non-polyhedral meshes creted in Gambit/Gmsh)
  if(PROGRAM_NAME .eq. 'Generate' .or.  &
     PROGRAM_NAME .eq. 'Convert') then
    call Enlarge % Matrix_Int(Grid % cells_bnd_region, i=(/1,6/),  &
                                                       j=(/-nb_m,nc_m/))
  end if

  ! Array to hold boundary regions tags at boundary cells
  if(nb_m .gt. 0) then  ! this check is not crazy, a (sub)domain might
                        ! have no boundary cells, all could be periodic
                        ! or simply not touching any boundary regions
    call Enlarge % Array_Int(Grid % region % at_cell, i=(/-nb_m,-1/))
  end if

  ! Allocate cell-based porous regions
  call Enlarge % Array_Int(Grid % por, i=(/-nb_m,nc_m/))

  ! Allocate processor i.d.
  call Enlarge % Array_Int(Grid % Comm % cell_proc, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Int(Grid % Comm % cell_glo,  i=(/-nb_m,nc_m/))
  do c = -nb_m, nc_m
    Grid % Comm % cell_glo(c) = c
  end do

  ! Allocate thread i.d.
  call Enlarge % Array_Int(Grid % Omp % cell_thread, i=(/-nb_m,nc_m/))

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  if(PROGRAM_NAME(1:7) .ne. 'Process') then
    call Enlarge % Array_Int(Grid % new_c, i=(/-nb_m,nc_m/))
    call Enlarge % Array_Int(Grid % old_c, i=(/-nb_m,nc_m/))
  end if

  end subroutine
