!==============================================================================!
  subroutine Print_Generate_Statistics(Generate, Grid)
!------------------------------------------------------------------------------!
!>  Prints some statistical data about the Grid on the terminal.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod, only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Generate_Type) :: Generate  !! parent class
  type(Grid_Type)      :: Grid      !! grid being analyzed
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, numb, nonz, stencw
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Generate)
!==============================================================================!

  print 1, '#================================================'
  print 1, '# Generated grid statistics                      '
  print 1, '#------------------------------------------------'
  print 2, '# Number of nodes         :', Grid % n_nodes
  print 2, '# Number of cells         :', Grid % n_cells
  print 2, '# Number of faces         :', Grid % n_faces
  print 2, '# Number of shadow faces  :', Grid % n_shadows
  print 2, '# Number of boundary cells:', Grid % n_bnd_cells
  print 1, '#------------------------------------------------'
1 format (a50)
2 format (a28, i8)

  ! Find the number of non zero entries
  nonz = 0
  do i = 1,Grid % n_cells
    stencw = 1            ! it used to be zero
    do j = 1, 24
      if( Grid % cells_c(j,i) > 0 ) stencw = stencw + 1
    end do
    nonz = nonz + stencw
  end do

  print 3, '# Number of non zero matrix entries:',   nonz
  print 4, '# Average stencil size             :',   real(nonz)  &
                                                   / real(Grid % n_cells)
  print 5, '#------------------------------------------------'
3 format (a37, i8)
4 format (a37, f8.3)
5 format (a50)

  ! Neighbours
  do j = 1, 24
    numb = 0
    do i = 1, Grid % n_cells
      stencw=0
      do k = 1, 24
        if( Grid % cells_c(k,i)  > 0 ) stencw=stencw+1
      end do
      if(stencw .eq. j) numb=numb+1
    end do
    if(numb .ne. 0) then
      print 6, '# Number of cells with ', j, ' neighbours: ', numb
    end if
  end do

  ! Twins
  do j = 1, 8
    numb = 0
    do i = 1, Grid % n_nodes
      if(Grid % twin_n(i,0) .eq. j) numb=numb+1
    end do
    if(numb .ne. 0) then
      print 6, '# Number of nodes with ', j, ' twins     : ', numb
    end if 
  end do
6 format (a24, i3, a13, i8)

  print *, '#------------------------------------------------'

  end subroutine
