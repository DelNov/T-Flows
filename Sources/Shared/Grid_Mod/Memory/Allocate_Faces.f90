!==============================================================================!
  subroutine Allocate_Faces(Grid, nf, ns, margin)
!------------------------------------------------------------------------------!
!>  Allocates memory for face-based data (arrays and matrices), for geometrical
!>  (xf, yf, zf, dx, dy, dz, ...) and connectivity data (faces_c, faces_n, ...).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)              :: Grid    !! grid being expanded
  integer, intent(in)           :: nf      !! number of faces in the grid
  integer, intent(in)           :: ns      !! number of shadow faces
  integer, intent(in), optional :: margin  !! margin for allocation
!-----------------------------------[Locals]-----------------------------------!
  integer :: m
!==============================================================================!

  ! Generator has growing number of faces, don't set them here
  if(PROGRAM_NAME(1:6) .ne. 'Genera') then
    Grid % n_faces   = nf
    Grid % n_shadows = ns
  end if

  ! If face-based arrays are allocated and they
  ! are bigger than requested, get out of here
  if(allocated(Grid % faces_n_nodes)) then
    if(nf + ns .le. ubound(Grid % faces_n_nodes, 1)) then
      return
    end if
  end if

  ! Process the margin if specified
  m = 0
  if(present(margin)) then
    Assert(margin .ge. 0)
    m = m + margin
  end if

  ! I can't figure out how to print something meaningful in parallel here
  if(Sequential_Run()) then
    print '(a,2i9)', ' # Expanding memory for faces to size: ', nf, ns
  end if

  ! Number of nodes at each face (determines face's shape really)
  call Enlarge % Array_Int (Grid % faces_n_nodes, i=(/1,nf+ns+m/))

  ! Faces' nodes, neigboring cells and shadows
  call Enlarge % Matrix_Int(Grid % faces_n, i=(/1,3/), j=(/1,nf+ns+m/))
  call Enlarge % Matrix_Int(Grid % faces_c, i=(/1,2/), j=(/1,nf+ns+m/))
  call Enlarge % Array_Int (Grid % faces_s, i=(/1,nf+ns+m/))

  ! Face surface areas (si), total surface (s)
  ! and distances between cells (di)
  call Enlarge % Array_Real(Grid % sx, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % sy, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % sz, i=(/1,nf+ns+m/))

  call Enlarge % Array_Real(Grid % s, i=(/1,nf+ns+m/))

  call Enlarge % Array_Real(Grid % dx, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % dy, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % dz, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % d,  i=(/1,nf+ns+m/))

  ! Face center coordinates
  call Enlarge % Array_Real(Grid % xf, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % yf, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % zf, i=(/1,nf+ns+m/))

  ! Vectors connecting face center with face cell centers connection
  call Enlarge % Array_Real(Grid % rx, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % ry, i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % rz, i=(/1,nf+ns+m/))

  ! Fractional cell volumes around faces
  call Enlarge % Array_Real(Grid % dv1,  i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % dv2,  i=(/1,nf+ns+m/))

  ! Weight factors
  call Enlarge % Array_Real(Grid % f,  i=(/1,nf+ns+m/))
  call Enlarge % Array_Real(Grid % fw, i=(/1,nf+ns+m/))

  ! Allocate thread i.d.
  call Enlarge % Array_Int(Grid % Omp % face_thread, i=(/1,nf+ns+m/))

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  if(PROGRAM_NAME(1:6) .ne. 'Proces') then
    call Enlarge % Array_Int(Grid % new_f, i=(/1,nf+ns+m/))
    call Enlarge % Array_Int(Grid % old_f, i=(/1,nf+ns+m/))
  end if

  ! Array to hold boundary regions tags at faces
  call Enlarge % Array_Int(Grid % region % at_face, i=(/1,nf+ns+m/))

  ! Variable which used to be declared in Gen_Mod:
  if(PROGRAM_NAME(1:6) .eq. 'Genera') then
    call Enlarge % Matrix_Int(Grid % face_c_to_c, i=(/1,nf+ns+m/), j=(/1,2/))
  end if

  end subroutine
