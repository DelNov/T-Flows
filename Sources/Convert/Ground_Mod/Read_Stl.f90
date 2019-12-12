!==============================================================================!
  subroutine Ground_Mod_Read_Stl(ground)
!------------------------------------------------------------------------------!
!   Reads the STL file with definition of the ground                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Ground_Type) :: ground
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_in
  integer           :: fu, cnt, v
  real              :: x_max, x_min, y_max, y_min, z_max, z_min, t_area
!==============================================================================!

  print *, '#==============================================================='
  print *, '# Enter the name of the ground file (with ext.):'
  print *, '#---------------------------------------------------------------'
  read(*,*) name_in

  call File_Mod_Open_File_For_Reading(name_in, fu)

  !--------------------------!
  !   Count all the facets   !
  !--------------------------!
  rewind(fu)
  cnt = 0
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. 'endsolid') exit
    if(line % tokens(1) .eq. 'facet') cnt = cnt + 1
  end do
  print *, '# Number of facets defining the ground: ', cnt

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  ground % n_facets = cnt
  allocate(ground % facet(ground % n_facets))

  !---------------------------------!
  !   Read all vertex coordinates   !
  !---------------------------------!
  x_max  = -HUGE
  y_max  = -HUGE
  z_max  = -HUGE
  x_min  = +HUGE
  y_min  = +HUGE
  z_min  = +HUGE
  t_area = 0.0
  rewind(fu)
  cnt = 0
  do
    call File_Mod_Read_Line(fu)
    if(line % tokens(1) .eq. 'endsolid') exit
    if(line % tokens(1) .eq. 'facet') then
      cnt = cnt + 1
      call File_Mod_Read_Line(fu)              ! 'outer loop'
      do v = 1, 3
        call File_Mod_Read_Line(fu)            ! 'vertex ... 1'
        read(line % tokens(2), *) ground % facet(cnt) % x(v)
        read(line % tokens(3), *) ground % facet(cnt) % y(v)
        read(line % tokens(4), *) ground % facet(cnt) % z(v)
        x_max = max(x_max, ground % facet(cnt) % x(v))
        x_min = min(x_min, ground % facet(cnt) % x(v))
        y_max = max(y_max, ground % facet(cnt) % y(v))
        y_min = min(y_min, ground % facet(cnt) % y(v))
        z_max = max(z_max, ground % facet(cnt) % z(v))
        z_min = min(z_min, ground % facet(cnt) % z(v))
      end do
      call File_Mod_Read_Line(fu)              ! 'endloop'
      ground % facet(cnt) % area_z =                     &
        Math_Mod_Tri_Area_Z(ground % facet(cnt) % x(1),  &
                            ground % facet(cnt) % y(1),  &
                            ground % facet(cnt) % x(2),  &
                            ground % facet(cnt) % y(2),  &
                            ground % facet(cnt) % x(3),  &
                            ground % facet(cnt) % y(3))
      t_area = t_area + ground % facet(cnt) % area_z
    end if
  end do
  print *, '# Read all ground coordinates!'
  print '(a,6f8.4)', ' # Bounding box: ', x_min, y_min, z_min, x_max, y_max, z_max
  print '(a, f8.4)', ' # Total area:   ', t_area

  close(fu)

  end subroutine
