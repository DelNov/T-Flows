!==============================================================================!
  subroutine Pick_A_Test_Case(Pol,         &
                              icellgeom,   &
                              ifunc)
!------------------------------------------------------------------------------!
!   Defines test case for individual cell / polyhedron                         !
!
!   icellgeom =  1, cube
!             =  2, tetrahedron
!             =  3, dodecahedron
!             =  4, icosahedron
!             =  5, complex cell (18 faces and 32 vertices)
!             =101, distorted cube
!             =102, non-convex pentagonal pyramid
!             =103, non-convex cell gto by subtracting a pyramid from the cube
!             =104, stellated cube
!             =105, non-convex hexahedron
!             =106, stellated dodecahedron
!             =107, stellated icosahedron
!             =108, hollowed cube
!             =109, drilled cube
!             =110, zig-zag prism
!  ifunc      =  1, sphere with radious 0.325 centered at (0.5,0.5,0.5)
!             =  2, torus with major radius 0.2, minor radius 0.1 and
!                  centered at (0.5,0.5,0.5)
!             =  3, orthocircle surface centered at (1.25,1.25,1.25)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
  integer                :: icellgeom, ifunc
!-----------------------------------[Locals]-----------------------------------!
  integer :: iread_error, ip
  real    :: x, y, z
!==============================================================================!

  print *, 'icellgeom:', icellgeom
  print *, 'ifunc:',     ifunc
  print *, 'phiiso:',    Pol % phi_iso

  if(icellgeom .lt. 1 .or. icellgeom .gt. 110 .or. (icellgeom.gt.5.and.  &
       icellgeom.lt.101)) then
  print *,'-----------------------------------------------------'
  print *,'|---------------------------------------------------|'
  print *,'|*********** WARNING FOR CELL SELECTION ************|'
  print *,'|---------------------------------------------------|'
  print *,'-----------------------------------------------------'
     print *, "Choose an appropriate icellgeom value (",          &
                 "between 1 and 5 for convex cells or between ",  &
                 "101 and 110 for non-convex cells)."
     stop
  end if
  if(ifunc.lt.0.or.ifunc.gt.3) then
  print *,'-----------------------------------------------------'
  print *,'|---------------------------------------------------|'
  print *,'|**** WARNING FOR MATERIAL BODY SHAPE SELECTION ****|'
  print *,'|---------------------------------------------------|'
  print *,'-----------------------------------------------------'
     print *, "Choose an appropriate ifunc value (between ",  &
                 "0 and 3)."
     stop
  end if
  if(icellgeom .eq. 1) then
     call Pol % Create_Cube()
  else if(icellgeom .eq. 2) then
     call Pol % Create_Tetrahedron()
  else if(icellgeom .eq. 3) then
     call Pol % Create_Dodecahedron()
  else if(icellgeom .eq. 4) then
     call Pol % Create_Icosahedron()
  else if(icellgeom .eq. 5) then
     call Pol % Create_Complexcell()
  else if(icellgeom .eq. 101) then
     call Pol % Create_Distortedcube()
  else if(icellgeom .eq. 102) then
     call Pol % Create_Pentapyramid()
  else if(icellgeom .eq. 103) then
     call Pol % Create_Cutcube()
  else if(icellgeom .eq. 104) then
     call Pol % Create_Scube()
  else if(icellgeom .eq. 105) then
     call Pol % Create_Nchexahedron()
  else if(icellgeom .eq. 106) then
     call Pol % Create_Sdodecahedron()
  else if(icellgeom .eq. 107) then
     call Pol % Create_Sicosahedron()
  else if(icellgeom .eq. 108) then
     call Pol % Create_Hollowedcube()
  else if(icellgeom .eq. 109) then
     call Pol % Create_Drilledcube()
  else if(icellgeom .eq. 110) then
     call Pol % Create_Zigzagcell()
  end if

  !-------------------------!
  !   Define scalar field   !
  !-------------------------!
  if(ifunc .eq. 0) open(1,file='phi',status='old',action='read')

  do ip = 1, Pol % n_nodes
    x = Pol % nodes_xyz(ip, 1)
    y = Pol % nodes_xyz(ip, 2)
    z = Pol % nodes_xyz(ip, 3)
    if(ifunc .eq. 0) then
      read(1,*,iostat=iread_error) Pol % phi(ip)
      if(iread_error .ne. 0) Pol % phi(ip)=0.0
    else if(ifunc .eq. 1) then
      Pol % phi(ip) = Pol % Func_1(x, y, z)
    else if(ifunc .eq. 2) then
      Pol % phi(ip) = Pol % Func_2(x, y, z)
    else if(ifunc .eq. 3) then
      Pol % phi(ip) = Pol % Func_3(x, y, z)
    end if
  end do

  if(ifunc .eq. 0) close(1)

  end

