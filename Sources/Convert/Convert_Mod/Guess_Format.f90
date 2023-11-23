!==============================================================================!
  subroutine Guess_Format(Convert, file_name, file_format)
!------------------------------------------------------------------------------!
!>  Designed to determine the format of a grid file based on its contents.
!>  As a result, returns one of the following: 'UNKNOWN', 'FLUENT', 'GMSH' or
!>  'GAMBIT' or 'OBJ' as the grid file format.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type)        :: Convert      !! parent class
  character(SL), intent(in)  :: file_name    !! input file name
  character(SL), intent(out) :: file_format  !! output format of the file
!-----------------------------------[Locals]-----------------------------------!
  integer :: l, fu
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  !----------------------------!
  !   Initialize file_format   !
  !----------------------------!
  file_format = 'UNKNOWN'

  ! Open the file as it is
  call File % Open_For_Reading_Ascii(file_name, fu)

  !-----------------------------------------!
  !   Check if it is Fluent's file format   !
  !-----------------------------------------!
  if(file_format .eq. 'UNKNOWN') then
    rewind(fu)
    do l = 1, 64
      call File % Read_Line(fu)
      if(Line % tokens(1) .eq. '(10') then
        file_format = 'FLUENT'
        print *, '#================================='   // &
                 '=================================='
        print *, '# Based on the contents, you are'     // &
                 ' reading Fluent''s .msh file format'
        print *, '#---------------------------------'   // &
                 '----------------------------------'
        exit
      end if
    end do
  end if

  !---------------------------------------!
  !   Check if it is Gmsh's file format   !
  !---------------------------------------!
  if(file_format .eq. 'UNKNOWN') then
    rewind(fu)
    do l = 1, 64
      call File % Read_Line(fu)
      if(Line % tokens(1) .eq. '$MeshFormat') then
        file_format = 'GMSH'
        print *, '#================================='   // &
                 '=================================='
        print *, '# Based on the contents, you are'     // &
                 ' reading Gmsh''s .msh file format'
        print *, '#---------------------------------'   // &
                 '----------------------------------'
        exit
      end if
    end do
  end if

  !-----------------------------------------!
  !   Check if it is Gambit's file format   !
  !-----------------------------------------!
  if(file_format .eq. 'UNKNOWN') then
    rewind(fu)
    do l = 1, 64
      call File % Read_Line(fu)
      if(Line % tokens(1) .eq. 'NUMNP' .and.  &
         Line % tokens(2) .eq. 'NELEM' .and.  &
         Line % tokens(3) .eq. 'NGRPS') then
        file_format = 'GAMBIT'
        print *, '#================================='   // &
                 '=================================='
        print *, '# Based on the contents, you are'     // &
                 ' reading Gambit''s .neu file format'
        print *, '#---------------------------------'   // &
                 '----------------------------------'
        exit
      end if
    end do
  end if

  !-----------------------------------------!
  !   Check if it is Gambit's file format   !
  !-----------------------------------------!
  if(file_format .eq. 'UNKNOWN') then
    rewind(fu)

    ! Skip three lines ... there could be "matlib"
    ! line and there are at least three vertice
    do l = 1, 3
      call File % Read_Line(fu)
    end do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'v') then
      file_format = 'OBJ'
      print *, '#=================================='   // &
               '==================================='
      print *, '# Based on the contents, you are re'   // &
               'ading Wavefront''s .obj file format'
      print *, '#----------------------------------'   // &
               '-----------------------------------'
      goto 1  ! exit wouldn't work here
    end if
  end if
1 continue

  close(fu)

  !----------------------------!
  !   File format is unknown   !
  !----------------------------!
  if(file_format .eq. 'UNKNOWN') then
    print *, '#================================================'
    print *, '# ERROR: Unrecognized input file format; exiting!'
    print *, '#------------------------------------------------'
    stop
  end if

  end subroutine
