!==============================================================================!
  subroutine Distribute_Cell_Coords(Surf)
!------------------------------------------------------------------------------!
!>  This subroutine is designed for distributing surface mesh coordinates
!>  across all processors in a parallel computing environment using MPI. It
!>  ensures that each processor gets the necessary data for vertices of the
!>  surface mesh, facilitating calculations in distributed computing systems.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias setup:                                                             !
!     - Establishes aliases for the grid (Grid) and vertices (Vert), and the   !
!       number of vertices (nv). This simplification enhances code readability !
!       and navigation.                                                        !
!   * Sequential and parallel processing:                                      !
!     - Differentiates between sequential and parallel execution environments. !
!     - In a sequential run, directly assigns cell coordinates (x, y, z) from  !
!       the grid to each vertex based on its nearest cell.                     !
!     - In a parallel run, handles the distribution of cell coordinates across !
!       multiple processors.                                                   !
!   * Buffer initialization in parallel runs:                                  !
!     - Initializes buffer arrays (buff_x, buff_y, buff_z) to store cell       !
!       coordinates and a count array (buff_n) for each vertex across          !
!       different processors.                                                  !
!   * Accumulation of coordinates:                                             !
!     - Iterates through each vertex in a parallel environment.                !
!     - Accumulates the x, y, and z coordinates from the vertex's nearest      !
!       cell if the cell is present in the current processor.                  !
!     - Keeps track of the number of accumulations in buff_n.                  !
!   * Global summation:                                                        !
!     - Performs a global summation of the buffer arrays across all processors !
!       using Global % Sum_Real_Array and Global % Sum_Int_Array. This step    !
!       ensures that each vertex receives data aggregated from all processors. !
!   * Updating vertex coordinates:                                             !
!     - After the global summation, updates each vertex's cell coordinates     !
!       by averaging the accumulated data. This process is crucial for         !
!       maintaining consistent cell coordinate data across the mesh in a       !
!       parallel computing environment.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Vert_Type), pointer :: Vert(:)
  integer,         pointer :: nv
  integer                  :: v, c
!==============================================================================!

  ! Take aliases
  Grid => Surf % pnt_grid
  Vert => Surf % Vert
  nv   => Surf % n_verts

  if(Sequential_Run()) then
    do v = 1, nv
      c = Surf % Vert(v) % cell  ! get nearest cell
      Surf % Vert(v) % cell_x = Grid % xc(c)
      Surf % Vert(v) % cell_y = Grid % yc(c)
      Surf % Vert(v) % cell_z = Grid % zc(c)
    end do
  else
    do v = 1, nv

      Surf % buff_x(v) = 0.0
      Surf % buff_y(v) = 0.0
      Surf % buff_z(v) = 0.0
      Surf % buff_n(v) = 0

      c = Surf % Vert(v) % cell  ! get nearest cell
      if(c > 0) then             ! if cell is in this processor
        if(Cell_In_This_Proc(c)) then
          Surf % buff_x(v) = Grid % xc(c)
          Surf % buff_y(v) = Grid % yc(c)
          Surf % buff_z(v) = Grid % zc(c)
          Surf % buff_n(v) = Surf % buff_n(v) + 1
        end if
      end if

    end do

    call Global % Sum_Real_Array(nv, Surf % buff_x)
    call Global % Sum_Real_Array(nv, Surf % buff_y)
    call Global % Sum_Real_Array(nv, Surf % buff_z)
    call Global % Sum_Int_Array (nv, Surf % buff_n)

    do v = 1, nv
      Surf % Vert(v) % cell_x = Surf % buff_x(v) / Surf % buff_n(v)
      Surf % Vert(v) % cell_y = Surf % buff_y(v) / Surf % buff_n(v)
      Surf % Vert(v) % cell_z = Surf % buff_z(v) / Surf % buff_n(v)
    end do
  end if

  end subroutine
