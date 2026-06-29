!==============================================================================!
  subroutine Sort_Nodes_Boundary_First(Grid, n_bnd_nodes)
!------------------------------------------------------------------------------!
!>  Sort nodes by placing boundary nodes first.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type),  intent(inout) :: Grid  !! parent grid
  integer, optional, intent(out)   :: n_bnd_nodes
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, c2, i_nod, n, cnt_bnd_n
  integer, allocatable :: key(:)
!==============================================================================!

  Assert(PROGRAM_NAME .eq. 'Convert')

  Assert(allocated(Grid % xn))

  write(*,'(a)', advance='no') ' # Sorting nodes boundary first ...'


  ! Reserve memory for working array
  call Enlarge % Array_Int(key, (/1, Grid % n_nodes/))
  key(:) = HUGE_INT

  cnt_bnd_n = 0  ! initialize number of boundary nodes

  ! Browse through all faces
  do s = 1, Grid % n_faces
    c2 = Grid % faces_c(2,s)

    ! Check if this is a boundary face
    if(c2 < 0) then

      do i_nod = 1, Grid % faces_n_nodes(s)
        n = Grid % faces_n(i_nod, s)
        if(key(n) .eq. HUGE_INT) then
          cnt_bnd_n = cnt_bnd_n + 1
          key(n) = cnt_bnd_n
        end if
      end do
    end if
  end do

  ! The sorting can take place now
  call Grid % Sort_Nodes_By_Key(key)

  if(present(n_bnd_nodes)) then
    n_bnd_nodes = cnt_bnd_n
  end if

  print '(a)', ' done!'

  end subroutine
