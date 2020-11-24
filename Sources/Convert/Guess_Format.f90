!==============================================================================!
  subroutine Guess_Format(file_name, file_format)
!------------------------------------------------------------------------------!
!   Based on the file contents, returns one of the following:                  !
!   'UNKNOWN', 'GAMBIT_NEU' or 'GMSH_MSH' as the grid file format.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(in)  :: file_name
  character(SL), intent(out) :: file_format
!-----------------------------------[Locals]-----------------------------------!
  integer :: l, fu
!==============================================================================!

  !----------------------------!
  !   Initialize file_format   !
  !----------------------------!
  file_format = 'UNKNOWN'

  ! Open the file as it is
  call File_Mod_Open_File_For_Reading(file_name, fu)

  !-----------------------------------------!
  !   Check if it is Fluent's file format   !
  !-----------------------------------------!
  if(file_format .eq. 'UNKNOWN') then
    rewind(fu)
    do l = 1, 64
      call File_Mod_Read_Line(fu)
      if(line % tokens(1) .eq. '(10') then
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
      call File_Mod_Read_Line(fu)
      if(line % tokens(1) .eq. '$MeshFormat') then
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
      call File_Mod_Read_Line(fu)
      if(line % tokens(1) .eq. 'NUMNP' .and.  &
         line % tokens(2) .eq. 'NELEM' .and.  &
         line % tokens(3) .eq. 'NGRPS') then
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
