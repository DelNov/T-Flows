!==============================================================================!
  program Main_Gis_2_Gmsh
!------------------------------------------------------------------------------!
!   Program which converts buildings in GIS format (as provided by Mahir)      !
!   to GMSH format; which can be included in scripts for generating meshes.    !
!                                                                              !
!   Compile with:                                                              !
!                                                                              !
!   gfortran -o ../../Binaries/Gis_2_Gmsh ../Shared/File_Mod.f90 Gis_2_Gmsh.f90
!                                                                              !
!   When prompted for the file name without extension, just type "buildings".  !
!   The program will assume that input file (as provided by Mahir) has the     !
!   extension ".gis", and will create the file "buildings.gmsh".               !
!                                                                              !
!   This "buildings.gms" is included in the script "sarajevo_3.city.gmsh".     !
!                                                                              !
!------------------------------------------------------------------------------!
  use File_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: MAXP = 65536   ! maximum number of points
  integer, parameter :: MAXB = MAXP/4  ! maximum number of buildings
  real,    parameter :: SCAL = 1.0e-3  ! coordinate scale
  real,    parameter :: HUGE = 1.0e+32 ! reall big number
  integer            :: f_in, f_out    ! file units
  integer            :: n_line, n_build, id_s, id_e, b, n, n0, n1, n2
  real               :: x(MAXP), y(MAXP)        ! points' coordinates
  real               :: del_min(MAXP)           ! points' minimum delta
  integer            :: build_id    (MAXB),  &  ! building's i.d. (tag from GIS)
                        build_s_line(MAXB),  &  ! building's starting line
                        build_e_line(MAXB)      ! building's ending line
  real               :: build_height(MAXB)
  logical            :: end_reached
  real               :: height, del, dot, dx_21, dy_21, dx_01, dy_01, u, xn, yn
  character(len=80)  :: name_in, name_out
!==============================================================================!

  print *, '#============================================================='
  print *, '# Program: Gis_2_Gmsh                                         '
  print *, '#                                                             '
  print *, '# Purpose: translates file describing bulings from GIS to     '
  print *, '#          another file which can be included in a GMSH script'
  print *, '#                                                             '
  print *, '# It is assumed that the file from GIS has exension ".gis"    '
  print *, '# whereas the file created by the program has exension ".gmsh"'
  print *, '#-------------------------------------------------------------'
  print *, '# Enter the GIS file name without extension:'
  read(*,*) problem_name

  !--------------------------------------------------!
  !   Set input (GIS) and output (GMSH) file names   !
  !--------------------------------------------------!
  call File_Mod_Set_Name(name_in,  extension=".gis");
  call File_Mod_Set_Name(name_out, extension=".gmsh");

  !---------------------------------!
  !   Open input and output files   !
  !---------------------------------!
  call File_Mod_Open_File_For_Reading(name_in,  f_in)
  call File_Mod_Open_File_For_Writing(name_out, f_out)

  !-----------------------------------!
  !                                   !
  !   Phase I: read the input files   !
  !                                   !
  !-----------------------------------!

  !---------------------------------------------------!
  !   Read input file to deduce:                      !
  !   - number of lines in the file,                  !
  !   - number of buildings,                          !
  !   - start and end line for each building.         !
  !---------------------------------------------------!
  n_line  = 0
  n_build = 0
  do
    call File_Mod_Read_Line(f_in, end_reached)
    if(end_reached) goto 1

    n_line = n_line + 1

    ! Read i.d., x, y and height from current line
    read(line % tokens(1), *) id_e
    read(line % tokens(2), *) x(n_line)
    read(line % tokens(3), *) y(n_line)
    read(line % tokens(6), *) height

    ! First line; first building starts
    if(n_line .eq. 1) then
      id_s = id_e                     ! 1st line: start of the 1st building
      n_build = n_build + 1           ! increase building count
      build_id    (n_build) = id_s    ! store builing's i.d.
      build_height(n_build) = height  ! store builing's height
      build_s_line(n_build) = n_line  ! set start of the new building
    end if

    ! All lines except first: when i.d.'s change, new building is starting
    if(id_e .ne. id_s) then
      id_s = id_e                         ! set start for the new building
      n_build = n_build + 1               ! new buildings encountered
      build_id    (n_build)   = id_s      ! store builing's i.d.
      build_height(n_build)   = height    ! store builing's height
      build_s_line(n_build)   = n_line    ! set start of the new and ...
      build_e_line(n_build-1) = n_line-1  ! ... end of the previous building
    end if
  end do
1 continue
  build_e_line(n_build) = n_line  ! finish the last building up

  print *, '# Number of lines in input file:',     n_line
  print *, '# Number of buildings in input file:', n_build

  !-------------------------------------------------!
  !                                                 !
  !   Phase II: Find minimum distance from points   !
  !             to other points and buldings        !
  !                                                 !
  !-------------------------------------------------!

  !----------------------------------------------------------!
  !   For each point find minumum distance to other points   !
  !----------------------------------------------------------!
  del_min(:) = HUGE
  do n1 = 1, n_line
    do n2 = n1 + 1, n_line
      del = sqrt( (x(n1) - x(n2))**2 + (y(n1)-y(n2))**2 )
      del_min(n1) = min(del_min(n1), del)
      del_min(n2) = min(del_min(n2), del)
    end do
  end do

  !------------------------------------------------------------!
  !   For each point find minumum distance to other segments   !
  !------------------------------------------------------------!
  do n0 = 1, n_line
    do b = 1, n_build

      ! Range for n if from first to last point in the building
      do n = build_s_line(b), build_e_line(b)

        ! Define segment in the building
        if(n < build_e_line(b)) then
          n1 = n
          n2 = n + 1
        else
          n1 = build_e_line(b)
          n2 = build_s_line(b)
        end if

        ! Treat only segments which don't contain n0
        if(n0 .ne. n1 .and. n0 .ne. n2) then
          dx_01 = x(n0) - x(n1)
          dy_01 = y(n0) - y(n1)
          dx_21 = x(n2) - x(n1)
          dy_21 = y(n2) - y(n1)
          u = (dx_01 * dx_21 + dy_01 * dy_21) / (dx_21**2 + dy_21**2)
          if(u >= 0.0 .and. u <= 1.0) then
            xn = x(n1) + u * dx_21
            yn = y(n1) + u * dy_21
            del = sqrt( (xn-x(n0))**2 + (yn-y(n0))**2 )
            del_min(n0) = min(del_min(n0), del)
            del_min(n1) = min(del_min(n1), del)
            del_min(n2) = min(del_min(n2), del)
          end if
        end if
      end do
    end do
  end do

  !---------------------------------------!
  !                                       !
  !   Phase III: write to output script   !
  !                                       !
  !---------------------------------------!
  write(f_out, '(a)') '//--------------------------------------------'
  write(f_out, '(a)') '// Definition of buildings begins             '
  write(f_out, '(a)') '//--------------------------------------------'
  write(f_out, '(a14,i4,a1)')  'n_buildings = ', n_build, ';'
  write(f_out, *)     ! an empty line

  n_line = 0          ! reset line counter
  do b = 1, n_build

    ! Building's header
    write(f_out, '(a21,i6)')     '// Building with ID: ', build_id(b)
    write(f_out, '(a12,i4,a1)')  'b         = ', b, ';'
    write(f_out, '(a12,i4,a1)')  'node_b(b) = ', 1 + build_e_line(b)  &
                                                   - build_s_line(b), ';'
    ! Building's node coordinates
    do n = 1, 1 + build_e_line(b) - build_s_line(b)

      n_line = n_line + 1  ! increase the line counter
      write(f_out, '(a10,i2,a4,f8.5,a3,                        &
                     a10,i2,a4,f8.5,a3,                        &
                     a10,i2,a4,f8.5,a1)')                      &
        'x(MAXN*b +', n, ') = ', x      (n_line)*SCAL, ';  ',  &
        'y(MAXN*b +', n, ') = ', y      (n_line)*SCAL, ';  ',  &
        'd(MAXN*b +', n, ') = ', del_min(n_line)*SCAL, ';'

    end do

    ! Building's height
    write(f_out, '(a14,i4,a1)') 'height_b(b) = ', nint(build_height(b)), ';'
    write(f_out, *)     ! an empty line
  end do

  end
