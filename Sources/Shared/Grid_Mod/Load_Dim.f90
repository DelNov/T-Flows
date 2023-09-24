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
  integer       :: fu, real_prec, version
  character(SL) :: name_in, str1, str2
  integer       :: nc, nb, nf, nn, ns
!==============================================================================!

  call Profiler % Start('Load_Dim')

  !----------------------------!
  !     Read the file with     !
  !   geometrical dimensions   !
  !----------------------------!
  call File % Set_Name(name_in, processor=this_proc, extension='.dim',  &
                       domain=domain)
  call File % Open_For_Reading_Binary(name_in, fu, this_proc)

  !----------------------------------------------!
  !   Store rank (domain number) for this grid   !
  !----------------------------------------------!
  Grid % rank = 0
  if(present(domain)) Grid % rank = domain

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
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells
  nn = Grid % n_nodes
  nf = Grid % n_faces
  ns = Grid % n_shadows

  call File % Buffered_Read_Real_Array(fu, Grid % xn(1:nn))
  call File % Buffered_Read_Real_Array(fu, Grid % yn(1:nn))
  call File % Buffered_Read_Real_Array(fu, Grid % zn(1:nn))

  call File % Buffered_Read_Real_Array(fu, Grid % xc(-nb:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % yc(-nb:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % zc(-nb:nc))

  ! Why on earth do I save wall distance for boundary cells?
  call File % Buffered_Read_Real_Array(fu, Grid % wall_dist(-nb:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % vol(1:nc))

  call File % Buffered_Read_Real_Array(fu, Grid % ixx(1:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % iyy(1:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % izz(1:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % ixy(1:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % ixz(1:nc))
  call File % Buffered_Read_Real_Array(fu, Grid % iyz(1:nc))

  call File % Buffered_Read_Real_Array(fu, Grid % sx(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % sy(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % sz(1:nf+ns))

  call File % Buffered_Read_Real_Array(fu, Grid % dx(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % dy(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % dz(1:nf+ns))

  call File % Buffered_Read_Real_Array(fu, Grid % f(1:nf+ns))

  call File % Buffered_Read_Real_Array(fu, Grid % xf(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % yf(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % zf(1:nf+ns))

  call File % Buffered_Read_Real_Array(fu, Grid % rx(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % ry(1:nf+ns))
  call File % Buffered_Read_Real_Array(fu, Grid % rz(1:nf+ns))

  read(fu) Grid % per_x
  read(fu) Grid % per_y
  read(fu) Grid % per_z

  close(fu)

  call Profiler % Stop('Load_Dim')

  end subroutine
