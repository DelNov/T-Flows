!==============================================================================!
  subroutine Save_Dim(Grid, sub)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for writing a .dim file containing
!>  geometrical data ("dimensions") of a computational grid. This file format
!>  to T-Flows and stands for "cells, faces, and nodes.".  Depending on the
!>  values in the Grid's fields new_n, old_c, and old_f, it will write a whole
!>  domain (from Generate and Convert), or an individual a sub-domain (from
!>  Divide).  The .dim and .cfn files (described in functions dealing with .cfn
!>  files) constitute T-Flows' native file format.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * File Creation: Initiates by creating a .dim file designated for storing  !
!     grid dimensions.                                                         !
!   * Data Writing:                                                            !
!     - Saves essential grid data such as real precision and file version.     !
!     - Writes node coordinates for all nodes in the grid.                     !
!     - Saves cell center coordinates, including for boundary cells and the    !
!       main computational cells.                                              !
!     - Records wall distance and cell volumes, crucial for various            !
!       computational fluid dynamics calculations.                             !
!     - Saves cell inertia tensors.                                            !
!     - For each face in the grid, writes a comprehensive set of data,         !
!       including face center coordinates, directional vectors, areas, and     !
!       other related geometric properties.                                    !
!     - Handles the directionality of faces correctly to ensure accurate       !
!       representation of the grid geometry.                                   !
!     - Records periodicity information of the grid if applicable.             !
!   * Efficient Data Handling: Utilizes buffer arrays for efficient handling   !
!     of large sets of geometric data.                                         !
!   * Subdomain Considerations: Handles subdomain information for parallel     !
!     processing, ensuring that only relevant data for the current subdomain   !
!     is processed and saved.                                                  !
!   * Finalization: Closes the file after writing all necessary data,          !
!     completing the operation.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid      !! grid being processed
  integer, intent(in) :: sub(1:2)  !! subdomain and total number of subdomains
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, s, n, c1, c2, var, fu, i
  character(SL)     :: name_out
  real, allocatable :: buffer(:)
!==============================================================================!

  call Profiler % Start('Save_Dim')

  !--------------------------!
  !   Allocate local array   !
  !--------------------------!
  allocate(buffer(max( Grid % n_nodes * 3,                            &
                      (Grid % n_bnd_cells + Grid % n_cells + 1) * 3,  &
                       Grid % n_faces + Grid % n_shadows)))

  !----------------------!
  !   Create .dim file   !
  !----------------------!
  call File % Set_Name(name_out, processor=sub, extension='.dim')
  call File % Open_For_Writing_Binary(name_out, fu)

  !-------------------------!
  !   Save real precision   !
  !-------------------------!
  write(fu) RP

  !------------------------------!
  !   Save version of the file   !
  !------------------------------!
  write(fu) VERSION_DIM

  !----------------------!
  !   Node coordinates   !
  !----------------------!
  i = 0  ! you allocated enough memory for all three components, check up
  do var = 1, 3
    do n = 1, Grid % n_nodes
      if(Grid % new_n(n) > 0) then
        i=i+1
        if(var .eq. 1) buffer(i) = Grid % xn(n)
        if(var .eq. 2) buffer(i) = Grid % yn(n)
        if(var .eq. 3) buffer(i) = Grid % zn(n)
      end if
    end do
  end do
  write(fu) buffer(1:i)

  !-----------------------------!
  !   Cell center coordinates   !
  !-----------------------------!
  i = 0  ! you allocated enough memory for all three components, check up
  do var = 1, 3
    do c = -Grid % n_bnd_cells, Grid % n_cells
      if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
        i=i+1
        if(var .eq. 1) buffer(i) = Grid % xc(Grid % old_c(c))
        if(var .eq. 2) buffer(i) = Grid % yc(Grid % old_c(c))
        if(var .eq. 3) buffer(i) = Grid % zc(Grid % old_c(c))
      end if
    end do
  end do
  write(fu) buffer(1:i)

  !-------------------!
  !   Wall distance   !  (buffered together with the next)
  !-------------------!
  i = 0  ! you allocated enough memory for all three components, check up
  do c = -Grid % n_bnd_cells, Grid % n_cells
    if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
      i=i+1;  buffer(i) = Grid % wall_dist(Grid % old_c(c))
    end if
  end do

  !------------------!
  !   Cell volumes   !  (buffered together with the previous)
  !------------------!
  do c = 1, Grid % n_cells
    if(Grid % old_c(c) .ne. 0) then
      i=i+1;  buffer(i) = Grid % vol(Grid % old_c(c))
    end if
  end do
  write(fu) buffer(1:i)

  !--------------------------!
  !   Cell inertia tensors   !
  !--------------------------!
  do var = 1, 6
    i = 0
    do c = 1, Grid % n_cells
      if(Grid % old_c(c) .ne. 0 .or. c .eq. 0) then
        i=i+1
        if(var .eq. 1) buffer(i) = Grid % ixx(Grid % old_c(c))
        if(var .eq. 2) buffer(i) = Grid % iyy(Grid % old_c(c))
        if(var .eq. 3) buffer(i) = Grid % izz(Grid % old_c(c))
        if(var .eq. 4) buffer(i) = Grid % ixy(Grid % old_c(c))
        if(var .eq. 5) buffer(i) = Grid % ixz(Grid % old_c(c))
        if(var .eq. 6) buffer(i) = Grid % iyz(Grid % old_c(c))
      end if
    end do
    write(fu) buffer(1:i)
  end do

  !-----------!
  !   Faces   !
  !-----------!
  do var = 1, 13
    i = 0
    do s = 1, Grid % n_faces + Grid % n_shadows
      if(Grid % old_f(s) .ne. 0) then
        c1 = Grid % faces_c(1, Grid % old_f(s))
        c2 = Grid % faces_c(2, Grid % old_f(s))
        if(Grid % new_c(c2) < 0 .or. Grid % new_c(c1) < Grid % new_c(c2)) then
          i=i+1
          if(var .eq.  1)  buffer(i) = Grid % sx(Grid % old_f(s))
          if(var .eq.  2)  buffer(i) = Grid % sy(Grid % old_f(s))
          if(var .eq.  3)  buffer(i) = Grid % sz(Grid % old_f(s))
          if(var .eq.  4)  buffer(i) = Grid % dx(Grid % old_f(s))
          if(var .eq.  5)  buffer(i) = Grid % dy(Grid % old_f(s))
          if(var .eq.  6)  buffer(i) = Grid % dz(Grid % old_f(s))
          if(var .eq.  7)  buffer(i) = Grid % f (Grid % old_f(s))
          if(var .eq.  8)  buffer(i) = Grid % xf(Grid % old_f(s))
          if(var .eq.  9)  buffer(i) = Grid % yf(Grid % old_f(s))
          if(var .eq. 10)  buffer(i) = Grid % zf(Grid % old_f(s))
          if(var .eq. 11)  buffer(i) = Grid % rx(Grid % old_f(s))
          if(var .eq. 12)  buffer(i) = Grid % ry(Grid % old_f(s))
          if(var .eq. 13)  buffer(i) = Grid % rz(Grid % old_f(s))
        else
          i=i+1
          if(var .eq.  1)  buffer(i) = -Grid % sx(Grid % old_f(s))
          if(var .eq.  2)  buffer(i) = -Grid % sy(Grid % old_f(s))
          if(var .eq.  3)  buffer(i) = -Grid % sz(Grid % old_f(s))
          if(var .eq.  4)  buffer(i) = -Grid % dx(Grid % old_f(s))
          if(var .eq.  5)  buffer(i) = -Grid % dy(Grid % old_f(s))
          if(var .eq.  6)  buffer(i) = -Grid % dz(Grid % old_f(s))
          if(var .eq.  7)  buffer(i) = 1.0 - Grid % f (Grid % old_f(s))
          if(var .eq.  8)  buffer(i) =  Grid % xf(Grid % old_f(s))
          if(var .eq.  9)  buffer(i) =  Grid % yf(Grid % old_f(s))
          if(var .eq. 10)  buffer(i) =  Grid % zf(Grid % old_f(s))
          if(var .eq. 11)  buffer(i) =  Grid % rx(Grid % old_f(s))
          if(var .eq. 12)  buffer(i) =  Grid % ry(Grid % old_f(s))
          if(var .eq. 13)  buffer(i) =  Grid % rz(Grid % old_f(s))
        end if
      end if
    end do
    write(fu) buffer(1:i)
  end do

  !-----------------!
  !   Periodicity   !
  !-----------------!
  write(fu) Grid % per_x
  write(fu) Grid % per_y
  write(fu) Grid % per_z

  close(fu)

  call Profiler % Stop('Save_Dim')

  end subroutine
