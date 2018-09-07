!==============================================================================!
  subroutine Load_Geo(grid, this_proc)
!------------------------------------------------------------------------------!
!   Reads:  name.geo                                                           !
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: this_proc
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, s
  character(len=80) :: name_in  
!==============================================================================!

  !----------------------------!
  !     Read the file with     !
  !   geometrical quantities   !
  !----------------------------!
  call Name_File(this_proc, name_in, '.geo') 
  open(9, file=name_in, form='unformatted')
  if(this_proc < 2) print *, '# Reading the file: ', name_in

  read(9) (grid % xn(n), n = 1, grid % n_nodes)
  read(9) (grid % yn(n), n = 1, grid % n_nodes) 
  read(9) (grid % zn(n), n = 1, grid % n_nodes)

  read(9) (grid % xc(c), c = 1, grid % n_cells)
  read(9) (grid % yc(c), c = 1, grid % n_cells) 
  read(9) (grid % zc(c), c = 1, grid % n_cells)

  read(9) (grid % xc(c), c=-1,-grid % n_bnd_cells,-1)  
  read(9) (grid % yc(c), c=-1,-grid % n_bnd_cells,-1)
  read(9) (grid % zc(c), c=-1,-grid % n_bnd_cells,-1) 

  read(9) (grid % vol  (c),     c = 1, grid % n_cells)
  read(9) (grid % delta(c),     c = 1, grid % n_cells)
  read(9) (grid % wall_dist(c), c = 1, grid % n_cells)

  read(9) (grid % sx(s), s = 1, grid % n_faces)
  read(9) (grid % sy(s), s = 1, grid % n_faces)
  read(9) (grid % sz(s), s = 1, grid % n_faces)

  read(9) (grid % dx(s), s = 1, grid % n_faces)
  read(9) (grid % dy(s), s = 1, grid % n_faces)
  read(9) (grid % dz(s), s = 1, grid % n_faces)

  read(9) (grid % f(s), s = 1, grid % n_faces)

  read(9) (grid % xf(s), s = 1, grid % n_faces)
  read(9) (grid % yf(s), s = 1, grid % n_faces)
  read(9) (grid % zf(s), s = 1, grid % n_faces)

  close(9) 

  end subroutine
