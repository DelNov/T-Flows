!==============================================================================!
  subroutine Grid_Mod_Load_Dim(grid, this_proc, domain)
!------------------------------------------------------------------------------!
!   Reads file with grid dimensions (.dim, used to be .geo)                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  integer           :: this_proc
  integer, optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, fu
  character(SL) :: name_in
!==============================================================================!

  !----------------------------!
  !                            !
  !     Read the file with     !
  !   geometrical dimensions   !
  !                            !
  !----------------------------!
  call File_Mod_Set_Name(name_in, processor=this_proc, extension='.dim',  &
                         domain=domain)
  call File_Mod_Open_File_For_Reading_Binary(name_in, fu, this_proc)

  read(fu) (grid % xn(n), n = 1, grid % n_nodes)
  read(fu) (grid % yn(n), n = 1, grid % n_nodes)
  read(fu) (grid % zn(n), n = 1, grid % n_nodes)

  read(fu) (grid % xc(c), c = -grid % n_bnd_cells, grid % n_cells)
  read(fu) (grid % yc(c), c = -grid % n_bnd_cells, grid % n_cells)
  read(fu) (grid % zc(c), c = -grid % n_bnd_cells, grid % n_cells)

  read(fu) (grid % wall_dist(c), c = -grid % n_bnd_cells, grid % n_cells)
  read(fu) (grid % vol(c), c = 1, grid % n_cells)

  read(fu) (grid % sx(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % sy(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % sz(s), s = 1, grid % n_faces + grid % n_shadows)

  read(fu) (grid % dx(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % dy(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % dz(s), s = 1, grid % n_faces + grid % n_shadows)

  read(fu) (grid % f(s), s = 1, grid % n_faces + grid % n_shadows)

  read(fu) (grid % xf(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % yf(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % zf(s), s = 1, grid % n_faces + grid % n_shadows)

  read(fu) (grid % rx(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % ry(s), s = 1, grid % n_faces + grid % n_shadows)
  read(fu) (grid % rz(s), s = 1, grid % n_faces + grid % n_shadows)

  read(fu) grid % per_x
  read(fu) grid % per_y
  read(fu) grid % per_z

  close(fu)

  end subroutine
