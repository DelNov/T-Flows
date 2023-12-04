!==============================================================================!
  subroutine Load_Cfn(Grid, procs, domain)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to read and load grid connectivity data from a
!>  .cfn file, essentially performing the opposite operation of the Save_Cfn
!>  subroutine.  This file format is specific to T-Flows and stands for "cells,
!>  faces, and nodes.".  Depending if it is called from Divide, Process in
!>  sequential mode or Process in parallel mode, it reads a whole domain
!>  or an individual a sub-domain.  The .cfn and .dim files (described in
!>  functions dealing with .dim files) constitute T-Flows' native file format.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * File reading setup: The subroutine begins by setting up the file name    !
!     for reading, considering processor information and subdomain ranking.    !
!   * Reading basic grid information:                                          !
!     - Reads the real precision used in the file and checks for compatibility !
!       with the current program's precision settings, halting if there's a    !
!       mismatch.                                                              !
!     - Reads the version of the .cfn file to ensure compatibility with the    !
!       current code version.                                                  !
!     - Retrieves the number of nodes, cells (including boundary cells), faces !
!       (including shadow faces), and boundary regions in the grid.            !
!   * Grid initialization and memory allocation: Allocates memory for nodes,   !
!     cells, faces, and boundary regions, preparing the Grid to store data.    !
!   * Reading detailed grid data:                                              !
!     - Reads the grid name and the list of boundary condition names.          !
!     - Retrieves global numbers for nodes.                                    !
!     - For cells, it reads the number of nodes and faces per cell, cells'     !
!       nodes, faces, processor IDs, global indices, and porosity regions.     !
!     - For faces, it reads the number of nodes per face, faces' nodes, cells  !
!       associated with each face, and shadows for each face. It performs      !
!       error checks to ensure data integrity.                                 !
!   * Boundary information:                                                    !
!     - Reads the boundary cell data, including information about the type of  !
!       boundary conditions applied.                                           !
!   * Error checking and validation: Throughout the process, it performs       !
!     several error checks to ensure the integrity of the read data. These     !
!     include checks for zero counts of nodes or faces for cells and faces,    !
!     zero indices for nodes in cells and faces, and proper linking between    !
!     real and shadow faces.                                                   !
!   * Finalizing: Once all data is read and validated, the file is closed,     !
!     and the subroutine concludes.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid        !! grid under consideration
  integer, intent(in) :: procs(1:2)  !! this and number of processors
  integer, optional   :: domain      !! rank of this subdomain
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, s, n, ss, sr, fu, real_prec, version
  character(SL)        :: name_in, str, str1, str2
  integer              :: nc, nb, nf, nn, ns, tot, i
  integer, allocatable :: i_buffer(:)
!==============================================================================!

  call Profiler % Start('Load_Cfn')

  !-------------------------------!
  !     Read the file with the    !
  !   connections between cells   !
  !-------------------------------!
  call File % Set_Name(name_in,           &
                       processor=procs,   &
                       extension='.cfn',  &
                       domain=domain)
  call File % Open_For_Reading_Binary(name_in, fu)

  !-------------------------!
  !   Read real precision   !
  !-------------------------!
  read(fu) real_prec

  if(real_prec .ne. RP) then
    if(RP .eq. 8 .and. real_prec .eq. 4) then
      call Message % Error(64,                                                 &
                       'Input files were saved in single precision, but '   // &
                       'this program is compiled in double. Decide in   '   // &
                       'which precision you want to work, recompile the '   // &
                       'programs (option REAL=double/single in make), '     // &
                       'convert or generate the meshes again and re-run.',     &
                       one_proc = .true.)
    else
      call Message % Error(64,                                                 &
                       'Input files were saved in double precision, but '   // &
                       'this program is compiled in single. Decide in   '   // &
                       'which precision you want to work, recompile the '   // &
                       'programs (option REAL=double/single in make), '     // &
                       'convert or generate the meshes again and re-run.',     &
                       one_proc = .true.)
    end if
  end if

  !------------------------------!
  !   Read version of the file   !
  !------------------------------!
  read(fu) version

  write(str1, '(i0.0)')  version
  write(str2, '(i0.0)')  VERSION_CFN

  ! Old version, but still backward compatible
  if(version .eq. 202304) then
    call Message % Warning(72,                                                 &
                 'You seem to be reading version '//trim(str1)//' of the  '//  &
                 '.cfn file, but the latest version is '//trim(str2)//'.  '//  &
                 'Porosity regions are not written in '//trim(str1)//'    '//  &
                 'and the code will set them all to zero.',                    &
                 one_proc = .true.)

  ! Old version, but incompatible
  else if(version .ne. VERSION_CFN) then
    call Message % Error(72,                                                   &
                 'You seem to be reading wrong version of the .cfn file.  '//  &
                 'The version you are reading is '//trim(str1)//' but the '//  &
                 'code expects version '//trim(str2)//'. Re-generate or   '//  &
                 'convert again the grids (and divide them if you run in  '//  &
                 'parallel).', one_proc = .true.)
  end if

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  read(fu) Grid % n_nodes
  read(fu) Grid % n_cells              ! number of cells including buffer
  read(fu) Grid % n_bnd_cells          ! number of boundary cells
  read(fu) Grid % n_faces              ! number of faces (with buffer faces)
  read(fu) Grid % n_shadows            ! number of shadow faces
  read(fu) Grid % n_bnd_regions        ! number of boundary conditions

  !-------------------------------------!
  !   Does grid have polyhedral cells   !
  !-------------------------------------!
  read(fu) Grid % polyhedral

  ! Allocate memory =--> carefull, there is no checking!
  call Grid % Allocate_Nodes  (Grid % n_nodes)
  call Grid % Allocate_Cells  (Grid % n_cells, Grid % n_bnd_cells)
  call Grid % Allocate_Faces  (Grid % n_faces, Grid % n_shadows)
  call Grid % Allocate_Regions(Grid % n_bnd_regions)

  ! Allocate buffer too, just for kicks
  allocate(i_buffer(Grid % n_faces + Grid % n_shadows))

  !-----------------!
  !   Domain name   !
  !-----------------!
  read(fu) Grid % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = Boundary_Regions()
    read(fu) Grid % region % name(n)
  end do

  ! Abbreviate from here on
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells
  nn = Grid % n_nodes
  nf = Grid % n_faces
  ns = Grid % n_shadows

  !--------------------------!
  !   Nodes global numbers   !
  !--------------------------!
  call File % Buffered_Read_Int_Array(fu, Grid % Comm % node_glo(1:nn))

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  call File % Buffered_Read_Int_Array(fu, Grid % cells_n_nodes(-nb:nc))

  ! Error trap for number of faces for each cell ...
  ! ... but also adjust array dimension properly ...
  ! ... and estimate the necessary buffer size (tot)
  tot = 0
  do c = -nb, nc
    if(c .ne. 0) then
      if(Grid % cells_n_nodes(c) .eq. 0) then
        write(str, '(i0.0)') c
        call Message % Error(72,                                           &
                   'Number of nodes is zero at cell: '//trim(str)//'. '//  &
                   'This is critical.  Exiting!',                          &
                   file=__FILE__, line=__LINE__)
      else
        call Enlarge % Matrix_Int(Grid % cells_n,  &
                            i=abs(Grid % cells_n_nodes(c)))
        tot = tot + abs(Grid % cells_n_nodes(c))
      end if
    end if
  end do
  call Enlarge % Array_Int(i_buffer, tot)

  ! Cells' nodes
  read(fu) i_buffer(1:tot)  ! guzzle the whole buffer at once
  i = 0
  do c = -nb, nc
    do n = 1, abs(Grid % cells_n_nodes(c))
      i=i+1;  Grid % cells_n(n, c) = i_buffer(i)
    end do
  end do

  ! Error trap for cells' nodes
  do c = -nb, nc
    do n = 1, abs(Grid % cells_n_nodes(c))
      if(Grid % cells_n(n, c) .eq. 0) then
        write(str, '(i0.0)') c
        call Message % Error(72,                                      &
                   'Node index is zero at cell: '//trim(str)//'. '//  &
                   'This error is critical.  Exiting!',               &
                   file=__FILE__, line=__LINE__)
      end if
    end do
  end do

  ! Number of faces for each cell
  call File % Buffered_Read_Int_Array(fu, Grid % cells_n_faces(-nb:nc))

  ! Error trap for number of faces for each cell ...
  ! ... but also adjust array dimension properly ...
  ! ... and estimate the necessary buffer size (tot)
  tot = 0
  do c = -nb, nc
    if(c .ne. 0) then
      if(Grid % cells_n_faces(c) .eq. 0) then
        write(str, '(i0.0)') c
        call Message % Error(72,                                           &
                   'Number of faces is zero at cell: '//trim(str)//'. '//  &
                   'This is critical.  Exiting!',                          &
                   file=__FILE__, line=__LINE__)
      else
        call Enlarge % Matrix_Int(Grid % cells_f,  &
                                i=Grid % cells_n_faces(c))
        tot = tot + Grid % cells_n_faces(c)
      end if
    end if
  end do
  call Enlarge % Array_Int(i_buffer, tot)

  ! Cells' faces
  read(fu) i_buffer(1:tot)  ! guzzle the whole buffer at once
  i = 0
  do c = -nb, nc
    do s = 1, Grid % cells_n_faces(c)
      i=i+1;  Grid % cells_f(s, c) = i_buffer(i)
    end do
  end do

  ! Error trap for cells' faces
  do c = -nb, nc
    do s = 1, Grid % cells_n_faces(c)
      if(Grid % cells_f(s, c) .eq. 0) then
        write(str, '(i0.0)') c
        call Message % Error(72,                                      &
                   'Face index is zero at cell: '//trim(str)//'. '//  &
                   'This error is critical.  Exiting!',               &
                   file=__FILE__, line=__LINE__)
      end if
    end do
  end do

  ! Cells' processor ids
  call File % Buffered_Read_Int_Array(fu, Grid % Comm % cell_proc(-nb:nc))

  ! Cells' global indices
  call File % Buffered_Read_Int_Array(fu, Grid % Comm % cell_glo(-nb:nc))

  ! Cells' porosity regions
  if(version .eq. VERSION_CFN) then
    call File % Buffered_Read_Int_Array(fu, Grid % por(-nb:nc))
  else if(version .eq. 202304) then
    Grid % por(-nb:nc) = 0
  end if

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  call File % Buffered_Read_Int_Array(fu, Grid % faces_n_nodes(1:ns+nf))

  ! Error trap for number of faces for each cell ...
  ! ... but also adjust array dimension properly ...
  ! ... and estimate the necessary buffer size (tot)
  tot = 0
  do s = 1, nf + ns
    if(Grid % faces_n_nodes(s) .eq. 0) then
      write(str, '(i0.0)') s
      call Message % Error(72,                                           &
                 'Number of nodes is zero at face: '//trim(str)//'. '//  &
                 'This is critical.  Exiting!',                          &
                 file=__FILE__, line=__LINE__)
    else
      call Enlarge % Matrix_Int(Grid % faces_n, i=Grid % faces_n_nodes(s))
      tot = tot + Grid % faces_n_nodes(s)
    end if
  end do
  call Enlarge % Array_Int(i_buffer, tot)

  ! Faces' nodes
  read(fu) i_buffer(1:tot)  ! guzzle the whole buffer at once
  i = 0
  do s = 1, nf + ns
    do n = 1, Grid % faces_n_nodes(s)
      i=i+1;  Grid % faces_n(n, s) = i_buffer(i)
    end do
  end do

  ! Error trap for faces' nodes
  do s = 1, nf + ns
    do n = 1, Grid % faces_n_nodes(s)
      if(Grid % faces_n(n, s) .eq. 0) then
        write(str, '(i0.0)') s
        call Message % Error(72,                                      &
                   'Node index is zero at face: '//trim(str)//'. '//  &
                   'This error is critical.  Exiting!',               &
                   file=__FILE__, line=__LINE__)
      end if
    end do
  end do

  ! Faces' cells
  tot = (nf + ns) * 2;
  call Enlarge % Array_Int(i_buffer, tot)
  read(fu) i_buffer(1:tot)  ! guzzle the whole buffer at once
  i = 0
  do s = 1, nf + ns
    i=i+1;  Grid % faces_c(1, s) = i_buffer(i)
    i=i+1;  Grid % faces_c(2, s) = i_buffer(i)
  end do

  ! Error trap for faces' cells
  do s = 1, nf + ns
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    ! Check only if least one cell is in this processor
    ! (Meaning it is not a face entirelly in the buffer)
    if(Grid % Comm % cell_proc(c1) .eq. procs(1) .or.  &
       Grid % Comm % cell_proc(c2) .eq. procs(1)) then
      if( .not. (c1.eq.0 .and. c2.eq.0) ) then
        if(Grid % faces_c(1, s) .eq. 0) then
          write(str,  '(i0.0)') s
          write(str1, '(i0.0)') c1;  write(str2, '(i0.0)') c2;
          call Message % Error(72,                                            &
                     'Cell one is zero at face: '//trim(str)//' '//           &
                     'surrounded by cells '//trim(str1)//' and '//trim(str2)  &
                     //'. \n This error is critical.  Exiting!',              &
                     file=__FILE__, line=__LINE__)
        end if
        if(Grid % faces_c(2, s) .eq. 0) then
          write(str,  '(i0.0)') s
          write(str1, '(i0.0)') c1;  write(str2, '(i0.0)') c2;
          call Message % Error(72,                                            &
                     'Cell two is zero at face: '//trim(str)//' '//           &
                     'surrounded by cells '//trim(str1)//' and '//trim(str2)  &
                     //'. \n This error is critical.  Exiting!',              &
                     file=__FILE__, line=__LINE__)
        end if
      end if
    end if
  end do

  ! Faces' shadows
  call File % Buffered_Read_Int_Array(fu, Grid % faces_s(1:nf+ns))

  ! Error trap for shadows
  do ss = nf + 1, nf + ns
    sr = Grid % faces_s(ss)  ! real face from shadow data
    if(sr .eq. 0) then
      call Message % Error(72,                               &
                 'Shadow face points to zero face. '   //    &
                 'This error is critical.  Exiting!',        &
                 file=__FILE__, line=__LINE__)
    end if
    if(Grid % faces_s(sr) .ne. ss) then
      call Message % Error(60,                                            &
                 'Real and shadow face do not point to one another. ' //  &
                 ' \n This error is critical.  Exiting!',                 &
                 file=__FILE__, line=__LINE__)
    end if
  end do

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells (and all the faces)
  ! (This opens the oportunity to store bounary condition info in ...
  !  ... the faces thus ridding us of the "if(c2 < 0) then" checks)
  allocate (Grid % region % at_cell(-nb:-1))
  call File % Buffered_Read_Int_Array(fu, Grid % region % at_cell(-nb:-1))

  allocate (Grid % region % at_face(1:nf))
  Grid % region % at_face(1:nf) = 0

  call Grid % Determine_Regions_Ranges()

  close(fu)

  call Profiler % Stop('Load_Cfn')

  end subroutine
