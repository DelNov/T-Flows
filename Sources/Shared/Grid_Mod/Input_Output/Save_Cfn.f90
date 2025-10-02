!==============================================================================!
  subroutine Save_Cfn(Grid,        &
                      sub,         &  ! subdomain
                      nn_sub,      &  ! number of nodes in the sub. 
                      nc_sub,      &  ! number of cells in the sub. 
                      nf_sub,      &  ! number of faces in the sub.
                      ns_sub,      &  ! number of shadow faces
                      nbc_sub)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for writing a .cfn file containing
!>  all connectivity data of a computational grid. This file format is specific
!>  to T-Flows and stands for "cells, faces, and nodes.".  Depending on the
!>  arguments sent to the function, it can write a whole domain (Generate and
!>  Convert), or an individual a sub-domain (from Divide).  The .cfn and .dim
!>  files (described in functions dealing with .dim files) constitute T-Flows'
!>  native file format.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Preparation and error checking: The subroutine starts by allocating      !
!     necessary memory and performing extensive error checking. It ensures the !
!     integrity of grid data, like the correctness of the number of nodes per  !
!     cell, faces per cell, and nodes per face. Any inconsistency in these     !
!     aspects triggers a critical error, halting the operation.                !
!   * File creation and data saving:                                           !
!     - It creates a .cfn file and begins by saving essential grid data,       !
!       including real precision, file version, and information about the grid !
!       such as the number of nodes, cells, faces, shadow faces, and boundary  !
!       cells.                                                                 !
!     - The presence of polyhedral cells in the grid is recorded.              !
!   * Writing grid and boundary condition names: The subroutine saves the grid !
!     name and a list of boundary condition names associated with the grid.    !
!   * Saving global numbers: It writes the global numbers of nodes, which are  !
!     important in parallel processing contexts to allow creation of patterns  !
!     for communication between sub-domains and consistent creation of backup  !
!     files.                                                                   !
!   * Cells, faces, and boundary data:                                         !
!     - For each cell, it saves the number of nodes and faces associated with  !
!       the cell, the global indices, processor IDs, and porosity regions.     !
!     - For each face, it stores the number of nodes, node indices, associated !
!       cells, and shadows. The order of nodes in a face is carefully handled  !
!       to maintain correct orientation.                                       !
!     - Physical boundary cell data are also saved, indicating the type of     !
!       boundary condition applied to each boundary cell.                      !
!   * Subdomain and parallel processing considerations:                        !
!     - The subroutine takes into account subdomain information (in parallel   !
!       processing scenarios), ensuring that only relevant data pertaining to  !
!       the current subdomain is processed and saved.                          !
!     - It uses efficient memory handling to cope with the potentially         !
!       large amount of data in a grid.                                        !
!   * Closing operations: Once all the necessary data is written, the file is  !
!     closed, and the subroutine concludes its operation.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid      !! grid being processed
  integer, intent(in) :: sub(1:2)  !! subdomain and total number of subdomains
  integer, intent(in) :: nn_sub    !! number of nodes in the (sub)domain
  integer, intent(in) :: nc_sub    !! number of cells in the (sub)domain
  integer, intent(in) :: nf_sub    !! number of faces in the (sub)domain
  integer, intent(in) :: ns_sub    !! number of shadows in the (sub)domain
  integer, intent(in) :: nbc_sub   !! number of boundary cells
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, n, i_nod, s, fu, c1, c2, ss, sr, max_n, i, j
  integer, allocatable :: faces_n(:), buffer(:)
  character(SL)        :: name_out, str, str1, str2
!==============================================================================!

  call Profiler % Start('Save_Cfn')

  ! Allocate memory for locals
  n = size(Grid % faces_n, 1)
  Assert(n .le. maxval(Grid % faces_n_nodes))
  allocate(faces_n(n))

  !--------------------------------------------------------!
  !                                                        !
  !   First all the error traps  (They are surprisingly    !
  !   fast compared to other sections of the subroutine)   !
  !                                                        !
  !--------------------------------------------------------!
  n = 0
  max_n = 0

  call Profiler % Start('Save_Cfn (error traps only)')

  ! Error trap for number of nodes for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      if(c .ne. 0) then
        if(Grid % cells_n_nodes(Grid % old_c(c)) .eq. 0) then
          write(str, '(i0.0)') Grid % old_c(c)
          call Message % Error(72,                                           &
                     'Number of nodes is zero at cell: '//trim(str)//'. '//  &
                     'This is critical.  Exiting!',                          &
                     file=__FILE__, line=__LINE__)
        end if
      end if
    end if
  end do

  ! Error trap for cells' nodes
  n = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      n = n + abs(Grid % cells_n_nodes(Grid % old_c(c)))
      do i_nod = 1, abs(Grid % cells_n_nodes(Grid % old_c(c)))
        if(Grid % new_n(Grid % cells_n(i_nod, Grid % old_c(c))) .eq. 0) then
          write(str, '(i0.0)') Grid % old_c(c)
          call Message % Error(72,                                      &
                     'Node index is zero at cell: '//trim(str)//'. '//  &
                     'This error is critical.  Exiting!',               &
                     file=__FILE__, line=__LINE__)
        end if
      end do
    end if
  end do
  max_n = max(max_n, n)

  ! Error trap for number of faces for each cell
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      if(c .ne. 0) then
        if(Grid % cells_n_faces(Grid % old_c(c)) .eq. 0) then
          write(str, '(i0.0)') Grid % old_c(c)
          call Message % Error(72,                                           &
                     'Number of faces is zero at cell: '//trim(str)//'. '//  &
                     'This is critical.  Exiting!',                          &
                     file=__FILE__, line=__LINE__)
        end if
      end if
    end if
  end do

  ! Error trap for cells' faces
  n = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      n = n + Grid % cells_n_faces(Grid % old_c(c))
      do s = 1, Grid % cells_n_faces(Grid % old_c(c))
        if(Grid % new_f(Grid % cells_f(s, Grid % old_c(c))) .eq. 0) then
          write(str, '(i0.0)') Grid % old_c(c)
          call Message % Error(72,                                      &
                     'Face index is zero at cell: '//trim(str)//'. '//  &
                     'This error is critical.  Exiting!',               &
                     file=__FILE__, line=__LINE__)
        end if
      end do
    end if
  end do
  max_n = max(max_n, n)

  ! Error trap for number of nodes for each face
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      if(Grid % faces_n_nodes(Grid % old_f(s)) .eq. 0) then
        write(str, '(i0.0)') Grid % old_f(s)
          call Message % Error(72,                                           &
                     'Number of nodes is zero at face: '//trim(str)//'. '//  &
                     'This is critical.  Exiting!',                          &
                     file=__FILE__, line=__LINE__)
      end if
    end if
  end do

  ! Error trap for faces' nodes
  n = 0
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      n = n + Grid % faces_n_nodes(Grid % old_f(s))
      do i_nod = 1, Grid % faces_n_nodes(Grid % old_f(s))
        if(Grid % new_n(Grid % faces_n(i_nod, Grid % old_f(s))) .eq. 0) then
        write(str, '(i0.0)') Grid % old_f(s)
          call Message % Error(72,                                      &
                     'Node index is zero at face: '//trim(str)//'. '//  &
                     'This error is critical.  Exiting!',               &
                     file=__FILE__, line=__LINE__)
        end if
      end do
    end if
  end do
  max_n = max(max_n, n)

  ! Error trap for faces' cells
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      c1 = Grid % faces_c(1, Grid % old_f(s))
      c2 = Grid % faces_c(2, Grid % old_f(s))

      ! Check only if least one cell is in this processor
      ! (Meaning it is not a face entirelly in the buffer)
      if(Grid % Comm % cell_proc(c1) .eq. sub(1) .or.  &
         Grid % Comm % cell_proc(c2) .eq. sub(1)) then
        if(Grid % new_c(c1) .eq. 0) then
          write(str,  '(i0.0)') Grid % old_f(s)
          write(str1, '(i0.0)') c1;  write(str2, '(i0.0)') c2;
          call Message % Error(72,                                            &
                     'Cell one is zero at face: '//trim(str)//' '//           &
                     'surrounded by cells '//trim(str1)//' and '//trim(str2)  &
                     //'. \n This error is critical.  Exiting!',              &
                     file=__FILE__, line=__LINE__)
        end if
        if(Grid % new_c(c2) .eq. 0) then
          write(str,  '(i0.0)') Grid % old_f(s)
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
  allocate(buffer(max_n))

  call Profiler % Stop('Save_Cfn (error traps only)')

  !-----------------------------------------------------!
  !                                                     !
  !   Enough trapping errors, create the file finally   !
  !                                                     !
  !-----------------------------------------------------!

  !----------------------!
  !   Create .cfn file   !
  !----------------------!
  call File % Set_Name(name_out,         &
                       processor=sub,    &
                       extension='.cfn')
  call File % Open_For_Writing_Binary(name_out, fu)

  !-------------------------!
  !   Save real precision   !
  !-------------------------!
  write(fu) RP

  !------------------------------!
  !   Save version of the file   !
  !------------------------------!
  write(fu) VERSION_CFN

  !-----------------------------------------------!
  !   Number of cells, boundary cells and faces   !
  !-----------------------------------------------!
  write(fu) nn_sub
  write(fu) nc_sub            ! new way: add buffer cells to cells
  write(fu) nbc_sub           ! number of boundary cells
  write(fu) nf_sub
  write(fu) ns_sub
  write(fu) Grid % n_bnd_regions  ! number of bounary conditions

  !-------------------------------------!
  !   Does grid have polyhedral cells   !
  !-------------------------------------!
  write(fu) Grid % polyhedral

  !------------------------!
  !   Domain (grid) name   !
  !------------------------!
  write(fu) Grid % name

  !------------------------------!
  !   Boundary conditions list   !
  !------------------------------!
  do n = Boundary_Regions()
    write(fu) Grid % region % name(n)
  end do

  !--------------------------!
  !   Nodes global numbers   !
  !--------------------------!
  i = 0
  do n = 1, Grid % n_nodes
    if(Grid % new_n(n) > 0) then
      i=i+1;  buffer(i) = Grid % Comm % node_glo(n)
    end if
  end do
  write(fu) buffer(1:i)

  !-----------!
  !   Cells   !  (including buffer cells)
  !-----------!

  ! Number of nodes for each cell
  i = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      i=i+1;  buffer(i) = Grid % cells_n_nodes(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  ! Cells' nodes
  i = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do j = 1, abs(Grid % cells_n_nodes(Grid % old_c(c)))
        i=i+1;  buffer(i) = Grid % new_n(Grid % cells_n(j, Grid % old_c(c)))
      end do
    end if
  end do
  write(fu) buffer(1:i)

  ! Number of faces for each cell
  i = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      i=i+1;  buffer(i) = Grid % cells_n_faces(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  ! Cells' faces.  They are still faces, kept in the same array.
  i = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      do s = 1, Grid % cells_n_faces(Grid % old_c(c))
        i=i+1;  buffer(i) = Grid % new_f(Grid % cells_f(s, Grid % old_c(c)))
      end do
    end if
  end do
  write(fu) buffer(1:i)

  ! Cells' processor ids
  i = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      i=i+1;  buffer(i) = Grid % Comm % cell_proc(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  ! Cells' global indices
  i = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      i=i+1;  buffer(i) = Grid % Comm % cell_glo(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  ! Cells' porosity regions
  i = 0
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      i=i+1;  buffer(i) = Grid % por(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  !-----------!
  !   Faces   !
  !-----------!

  ! Number of nodes for each face
  i = 0
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      i=i+1;  buffer(i) = Grid % faces_n_nodes(Grid % old_f(s))
    end if
  end do
  write(fu) buffer(1:i)

  ! Faces' nodes
  i = 0
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      c1 = Grid % faces_c(1, Grid % old_f(s))
      c2 = Grid % faces_c(2, Grid % old_f(s))

      ! Shorten the syntax a bit
      n = Grid % faces_n_nodes(Grid % old_f(s))

      ! Copy the nodes in local variable
      do i_nod = 1, n
        faces_n(i_nod) = Grid % new_n(Grid % faces_n(i_nod, Grid % old_f(s)))
      end do

      ! Write the local variable out taking care of the order
      if(Grid % new_c(c2) < 0 .or. Grid % new_c(c1) < Grid % new_c(c2)) then
        do j = 1, n  ! for buffers, I decided to use "j" and not "i_nod"
          i=i+1;  buffer(i) = faces_n(j)
        end do
      else
        do j = n, 1, -1
          i=i+1;  buffer(i) = faces_n(j)
        end do
      end if
    end if
  end do
  write(fu) buffer(1:i)

  ! Faces' cells
  i = 0
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(s) .ne. 0) then
      c1 = Grid % faces_c(1, Grid % old_f(s))
      c2 = Grid % faces_c(2, Grid % old_f(s))

      ! At least one cell is in this processor
      ! (Meaning it is not a face entirelly in the buffer)
      if(Grid % Comm % cell_proc(c1) .eq. sub(1) .or.  &
         Grid % Comm % cell_proc(c2) .eq. sub(1)) then
        if(Grid % new_c(c2) < 0 .or. Grid % new_c(c1) < Grid % new_c(c2)) then
          i=i+1;  buffer(i) = Grid % new_c(c1)
          i=i+1;  buffer(i) = Grid % new_c(c2)
        else
          i=i+1;  buffer(i) = Grid % new_c(c2)
          i=i+1;  buffer(i) = Grid % new_c(c1)
        end if

      ! Face is not active
      else
        i=i+1;  buffer(i) = 0
        i=i+1;  buffer(i) = 0
      end if
    end if
  end do
  write(fu) buffer(1:i)

  ! Faces' shadows
  i = 0
  do sr = 1, Grid % n_faces
    if(Grid % old_f(sr) .ne. 0) then         ! this face will be saved
      ss = Grid % faces_s(Grid % old_f(sr))  ! fetch its shadow

      if(ss .ne. 0) then  ! the saved face (sr) has a shadow (ss)
        i=i+1;  buffer(i) = Grid % new_f(ss)
      else
        i=i+1;  buffer(i) = 0
      end if
    end if
  end do

  ! Still shadows, don re-start the i
  do ss = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    if(Grid % old_f(ss) .ne. 0) then  ! this shadow face will be saved
      sr = Grid % faces_s(ss)         ! fetch its real counterpart

      if(sr .ne. 0) then  ! the saved face (sr) has a shadow (ss)
        i=i+1;  buffer(i) = Grid % new_f(sr)
      else
        call Message % Error(72,                               &
                   'Shadow face points to zero face. '   //    &
                   'This error is critical.  Exiting!',        &
                   file=__FILE__, line=__LINE__)
      end if
    end if
  end do
  write(fu) buffer(1:i)

  !--------------!
  !   Boundary   !
  !--------------!

  ! Physical boundary cells
  i = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % old_c(c) .ne. 0) then
      i=i+1;  buffer(i) = Grid % region % at_cell(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  close(fu)

  call Profiler % Stop('Save_Cfn')

  end subroutine
