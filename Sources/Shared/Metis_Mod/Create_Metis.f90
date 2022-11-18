!==============================================================================!
  subroutine Create_Metis(Metis, face_to_cell, n_parts)
!------------------------------------------------------------------------------!
!   Prepares all the arrays for a call to METIS function(s)                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Metis_Type)                                :: Metis
  integer, allocatable, intent(in), dimension(:,:) :: face_to_cell
  integer, optional,    intent(in)                 :: n_parts
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, i, v, v1, v2
!==============================================================================!

  Assert(n_parts > 0)

  !------------------------------------------------------------!
  !   Number of vertices and number of edges for first level   !
  !------------------------------------------------------------!
  Metis % n_verts = 0
  Metis % n_edges = 0
  do s = 1, size(face_to_cell, 2)
    if(face_to_cell(2, s) > 0) then
      Metis % n_edges = Metis % n_edges + 1
      Metis % n_verts = max(Metis % n_verts, face_to_cell(2, s))
    end if
  end do
  Assert(Metis % n_verts > 0)
  Assert(Metis % n_edges > 0)

  !---------------------------------------------------------------!
  !   Once n_verts(1) and n_edegs(1) are known, allocate memory   !
  !---------------------------------------------------------------!
  allocate(Metis % edges_v  ( 2,Metis % n_edges));  Metis % edges_v  (:,:) = 0
  allocate(Metis % edges_c  (   Metis % n_edges));  Metis % edges_c  (:)   = 0
  allocate(Metis % star_size(   Metis % n_verts));  Metis % star_size(:)   = 0
  allocate(Metis % star     (24,Metis % n_verts));  Metis % star     (:,:) = 0

  !----------------------------!
  !   Form edge connectivity   !
  !----------------------------!
  i = 0
  do s = 1, size(face_to_cell, 2)
    if(face_to_cell(2, s) > 0) then
      i = i + 1
      Metis % edges_v(1:2, i) = face_to_cell(1:2, s)
      Assert(Metis % edges_v(1,i) > 0)
      Assert(Metis % edges_v(2,i) > 0)
    end if
  end do

  !------------------------------!
  !   Form stars at this level   !
  !------------------------------!
  do s = 1, Metis % n_edges
    v1 = Metis % edges_v(1, s)
    v2 = Metis % edges_v(2, s)

    Metis % star_size(v1) = Metis % star_size(v1) + 1
    Metis % star_size(v2) = Metis % star_size(v2) + 1

    Metis % star(Metis % star_size(v1), v1) = v2
    Metis % star(Metis % star_size(v2), v2) = v1
  end do

  !---------------------------------------------------------------!
  !   Fill-up the structures needed to call METIS (row and col)   !
  !---------------------------------------------------------------!

  ! Fill up the rows
  allocate(Metis % row(Metis % n_verts + 1))
  Metis % row(1) = 0
  do v = 1, Metis % n_verts
    Metis % row(v+1) = Metis % row(v) + Metis % star_size(v)
  end do

  ! Fill up columns
  allocate(Metis % col(Metis % row(Metis % n_verts + 1)))
  do v = 1, Metis % n_verts
    do i = 1, Metis % star_size(v)
      ! In the line below -1 is used because METIS works from 0
      Metis % col(Metis % row(v) + i) = Metis % star(i, v) - 1
    end do
  end do

  !-------------------------------!
  !   Deal with constrains next   !
  !-------------------------------!
  Metis % n_constrains =  1
  allocate(Metis % imbalance(Metis % n_constrains))
  allocate(Metis % vert_weights(Metis % n_verts*Metis % n_constrains))
  allocate(Metis % edge_weights(Metis % row(Metis % n_verts+1)))
  allocate(Metis % vert_data(Metis % n_verts))
  allocate(Metis % part_weight(n_parts*Metis % n_constrains))
  Metis % imbalance(:)    = 1.001
  Metis % vert_weights(:) = 1
  Metis % edge_weights(:) = 1
  Metis % vert_data(:)    = 1
  Metis % part_weight(:)  = 1.0 / real(n_parts)

  end subroutine
