!==============================================================================!
  subroutine Find_Cells_Faces(Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is tasked with identifying all the faces that surround
!>  each cell in the computational grid, including boundary cells. It is
!>  comprehensive in that it also considers shadow faces
!------------------------------------------------------------------------------!
!   Process                                                                    !
!                                                                              !
!   * Initialization:                                                          !
!     - Resets the count of faces per cell to zero for all cells.              !
!   * Face Mapping:                                                            !
!     - Iteratively maps each face to its adjacent cells.                      !
!     - Includes logic to consider shadow faces alongside regular faces.       !
!     - Uses geometrical calculations to determine the proximity of shadow     !
!       faces to the cells, opting for the closer face when necessary.         !
!   * Finalization:                                                            !
!     - Concludes the mapping process and stops the profiling timer.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_cel, s, sh   ! counters
  real    :: dist_s, dist_sh
!==============================================================================!

  call Profiler % Start('Find_Cells_Faces')

  if(any(Grid % cells_n_faces(:) .ne. 0)) then
    print *, '# NOTE: Seems you are looking for cells'' faces'
    print *, '# although this information has already been found!'
    print *, '# No harm done, just a little note.'
  end if

  Grid % cells_n_faces(:) = 0

  do s = 1, Grid % n_faces
    do i_cel = 1, 2
      c = Grid % faces_c(i_cel, s)  ! would be c1 and c2 in most of the code

      Grid % cells_n_faces(c) = Grid % cells_n_faces(c) + 1
      call Enlarge % Matrix_Int(Grid % cells_f, i=Grid % cells_n_faces(c))
      Grid % cells_f(Grid % cells_n_faces(c), c) = s

      sh = Grid % faces_s(s)        ! take the shadow face
      if(sh > 0) then
        dist_s  = Math % Distance(                                  &
                       Grid % xc(c),  Grid % yc(c),  Grid % zc(c),  &
                       Grid % xf(s),  Grid % yf(s),  Grid % zf(s))
        dist_sh = Math % Distance(                                  &
                       Grid % xc(c),  Grid % yc(c),  Grid % zc(c),  &
                       Grid % xf(sh), Grid % yf(sh), Grid % zf(sh))
        if(dist_sh < dist_s) then
          Grid % cells_f(Grid % cells_n_faces(c), c) = s
        end if
      end if
    end do
  end do

  call Profiler % Stop('Find_Cells_Faces')

  end subroutine
