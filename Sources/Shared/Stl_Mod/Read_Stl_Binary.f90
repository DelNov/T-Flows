!==============================================================================!
  subroutine Read_Stl_Binary(Stl)
!------------------------------------------------------------------------------!
!>  This subroutine reads an STL (Stereolithography) file in binary format and
!>  populates an Stl_Type object with its contents. It parses facet normals and
!>  vertex coordinates to reconstruct the geometry defined in the STL file.
!------------------------------------------------------------------------------!
! implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type) :: Stl  !! parent Stl_Type object
!-----------------------------------[Locals]-----------------------------------!
  integer    :: fu, v, f, i_ver, n
  integer(1) :: byte
!==============================================================================!

  call File % Open_For_Reading_Binary(Stl % name, fu)

  ! First 80 characters are reserved and to be skiped
  do n = 1, 80
    read(fu) byte
  end do

  call File % Read_Binary_Int4_Array(fu, 1)
  f = int4_array(1)
  if(First_Proc()) print '(a,i9)', ' # Number of facets on the STL: ', f
  call Stl % Allocate_Stl(f)

  !------------------------------------------------!
  !   Read all facet normals, vertex coordinates   !
  !------------------------------------------------!
  do f = -Stl % n_bnd_cells, -1
    call File % Read_Binary_Real4_Array(fu, 3)
    Stl % nx(f) = real4_array(1)
    Stl % ny(f) = real4_array(2)
    Stl % nz(f) = real4_array(3)
    do i_ver = 1, 3
      v = (abs(f)-1) * 3 + i_ver   ! find global vertex number ...
      Stl % cells_n(i_ver, f) = v  ! ... and store it in cells_n
      call File % Read_Binary_Real4_Array(fu, 3)
      Stl % xn(v) = real4_array(1)
      Stl % yn(v) = real4_array(2)
      Stl % zn(v) = real4_array(3)
    end do
    read(fu) byte
    read(fu) byte
  end do
  if(First_Proc()) print '(a)', ' # Read all STL facets!'

  end subroutine
