!==============================================================================!
  subroutine Load_Dim(Grid, this_proc, domain)
!------------------------------------------------------------------------------!
!   Reads file with grid dimensions (.dim, used to be .geo)                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: this_proc
  integer, optional   :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, fu, real_prec, version
  character(SL) :: name_in, str1, str2
!==============================================================================!

  call Profiler % Start('Load_Dim')

  !----------------------------!
  !     Read the file with     !
  !   geometrical dimensions   !
  !----------------------------!
  call File % Set_Name(name_in, processor=this_proc, extension='.dim',  &
                       domain=domain)
  call File % Open_For_Reading_Binary(name_in, fu, this_proc)

  !-------------------------!
  !   Read real precision   !
  !-------------------------!
  read(fu) real_prec

  if(real_prec .ne. RP) then
    if(RP .eq. 8 .and. real_prec .eq. 4) then
      call Message % Error(64,                                                 &
                       'Input files were saved in single precision, but '   // &
                       'this program is compiled in double. Decide in   '   // &
                       'which precision you want to work, recompile the '   // &
                       'programs (option REAL=double/single in make), '     // &
                       'convert or generate the meshes again and re-run.',     &
                       one_proc = .true.)
    else
      call Message % Error(64,                                                 &
                       'Input files were saved in double precision, but '   // &
                       'this program is compiled in single. Decide in   '   // &
                       'which precision you want to work, recompile the '   // &
                       'programs (option REAL=double/single in make), '     // &
                       'convert or generate the meshes again and re-run.',     &
                       one_proc = .true.)
    end if
  end if

  !------------------------------!
  !   Read version of the file   !
  !------------------------------!
  read(fu) version

  if(version .ne. VERSION_DIM) then
    write(str1, '(i0.0)')  version
    write(str2, '(i0.0)')  VERSION_DIM
    call Message % Error(72,                                                   &
                 'You seem to be reading wrong version of the .dim file.  '//  &
                 'The version you are reading is '//trim(str1)//' but the '//  &
                 'code expects version '//trim(str2)//'. Re-generate or   '//  &
                 'convert again the grids (and divide them if you run in  '//  &
                 'parallel).', one_proc = .true.)
  end if

  !-------------------------!
  !   Read everything else  !
  !-------------------------!
  read(fu) (Grid % xn(n), n = 1, Grid % n_nodes)
  read(fu) (Grid % yn(n), n = 1, Grid % n_nodes)
  read(fu) (Grid % zn(n), n = 1, Grid % n_nodes)

  read(fu) (Grid % xc(c), c = -Grid % n_bnd_cells, Grid % n_cells)
  read(fu) (Grid % yc(c), c = -Grid % n_bnd_cells, Grid % n_cells)
  read(fu) (Grid % zc(c), c = -Grid % n_bnd_cells, Grid % n_cells)

  read(fu) (Grid % wall_dist(c), c = -Grid % n_bnd_cells, Grid % n_cells)
  read(fu) (Grid % vol(c), c = 1, Grid % n_cells)

  read(fu) (Grid % ixx(c), c = 1, Grid % n_cells)
  read(fu) (Grid % iyy(c), c = 1, Grid % n_cells)
  read(fu) (Grid % izz(c), c = 1, Grid % n_cells)
  read(fu) (Grid % ixy(c), c = 1, Grid % n_cells)
  read(fu) (Grid % ixz(c), c = 1, Grid % n_cells)
  read(fu) (Grid % iyz(c), c = 1, Grid % n_cells)

  read(fu) (Grid % sx(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % sy(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % sz(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % dx(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % dy(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % dz(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % f(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % xf(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % yf(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % zf(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % rx(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % ry(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % rz(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) Grid % per_x
  read(fu) Grid % per_y
  read(fu) Grid % per_z

  close(fu)

  call Profiler % Stop('Load_Dim')

  end subroutine
