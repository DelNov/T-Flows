!==============================================================================!
  subroutine Create_From_File(Stl, file_name)
!------------------------------------------------------------------------------!
!>  This subroutine constructs an STL object from a given file. It reads the
!>  contents of the STL file specified by 'file_name' and populates the Stl
!>  object with its data. The routine handles both binary and ASCII STL formats.
!>  It also merges duplicate vertices to streamline the object's structure.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)          :: Stl        !! parent Stl_Type object
  character(*), intent(in) :: file_name  !! file name
!-----------------------------------[Locals]-----------------------------------!
  integer :: v, f, i_ver, n
  integer :: body = 0
  logical :: und_1, und_2
!==============================================================================!

  ! Set the STL's name and open the STL file
  Stl % name = file_name

  !-----------------------!
  !   Read the STL file   !
  !-----------------------!
  if(File % Is_In_Binary_Format(Stl % name)) then
    call Stl % Read_Stl_Binary()
  else
    call Stl % Read_Stl_Ascii()

  end if

  !---------------------------------------------!
  !   This is brave: merge duplicate vertices   !
  !---------------------------------------------!
  allocate(Stl % new_n(Stl % n_nodes));  Stl % new_n(:) = 0
  allocate(Stl % old_n(Stl % n_nodes));  Stl % old_n(:) = 0
  call Stl % Merge_Duplicate_Nodes()

  !---------------------------------------------!
  !                                             !
  !---------------------------------------------!
  body = 0
  allocate(Stl % body_c(-Stl % n_bnd_cells:-1));  Stl % body_c(:) = 0.0
  allocate(Stl % body_n( Stl % n_nodes));         Stl % body_n(:) = 0.0

  do n = 1, 9999

    ! Find first und_1 facet
    und_1 = .false.
    do f = -Stl % n_bnd_cells, -1
      if(Stl % body_c(f) .eq. 0) then
        und_1 = .true.
        body = body + 1
        exit
      end if
    end do

    if(und_1) then

      ! Mark first facet and its vertices
      Stl % body_c(f) = body
      do i_ver = 1, 3
        v = Stl % cells_n(i_ver, f)
        Stl % body_n(v) = body
      end do

      ! Flood fill
1     continue
      und_2 = .false.
      do f = -Stl % n_bnd_cells, -1
        if(Stl % body_c(f) .eq. 0) then
          if(any(Stl % body_n(Stl % cells_n(1:3, f)) .eq. body)) then
            und_2 = .true.
            Stl % body_c(f) = body
            do i_ver = 1, 3
              v = Stl % cells_n(i_ver, f)
              Stl % body_n(v) = body
            end do
          end if
        end if
      end do
      if(und_2) goto 1
    else  ! no more undefiend facets
      exit
    end if

  end do

  Stl % n_boddies = maxval(Stl % body_c)
  if(First_Proc()) then
    print '(a,i3)', ' # Number of boddies in the STL file: ', Stl % n_boddies
  end if

  call Stl % Save_Debug_Vtu(append="body",                   &
                            scalar_name="body",              &
                            scalar_cell=real(Stl % body_c),  &
                            scalar_node=real(Stl % body_n),  &
                            plot_inside=.false.)


  end subroutine

