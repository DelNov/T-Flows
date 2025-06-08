!==============================================================================!
  subroutine Check_Matrix_Properties(Amg,         &
                                     a, ia, ja,   &  ! defining system
                                     icg, ifg)
!------------------------------------------------------------------------------!
!   Checks for several properties of a, ia, ja. in particular,
!   checks for symmetric storage of given matrix (i.e. l(i,j) is
!   stored in a iff l(j,i) is stored in a).
!
!   if storage is symmetric, pairs of zeroes are removed (if there
!   are any) and program execution con tinues.
!
!   If, however, storage of a is not symmetric, it is symmetrized.
!   if isym.eq.1, the missing transpose connections are copied.
!   if isym.eq.2, the storage of the matrix a is symmetrized (by
!   adding certain zero elements). then pairs of zeroes are removed
!   (if there are any) and program execution con tinues.
!
!   Arrays used for temporary storage: icg, ifg
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  real            :: a(:)
  integer         :: ia(:), ja(:)
  integer         :: icg(:), ifg(:)
!-----------------------------------[locals]-----------------------------------!
  real    :: anormm, anormp, asym, d, deps, rowsum
  integer :: i, i1, ishift, j, j1, j2, jnew
  integer :: napos, naoff, naneg, nazer, new, nna, nnu
  integer :: nda, ndicg
  logical :: found
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  nda   = size(a, 1)
  ndicg = size(icg, 1)

  !---------------------------------------------!
  !   Check if pointers ia, ja are reasonable   !
  !---------------------------------------------!
  if(ia(1) .ne. 1) then
    write(6, '(a)')  &
      ' *** error in check: pointer ia erroneous ***'
    Amg % ierr = AMG_ERR_IA_POINTER_INVALID
    return
  end if
  nnu = Amg % imax(1)
  if(nnu .ge. ndicg) then
    write(6, '(a)')  &
      ' *** error in check: ndw too small ***'
    Amg % ierr = AMG_ERR_DIM_ICG_TOO_SMALL
    return
  end if
  do i = 1, nnu
    icg(i) = 0
  end do

  !-----------------------------!
  !   Browse through unknowns   !
  !-----------------------------!
  do i = 1, nnu
    j1 = ia(i)
    j2 = ia(i+1)-1
    if (j2.lt.j1 .or. j2.gt.nda) then
      write(*,*) 'error, (j2.lt.j1 .or. j2.gt.nda) '
      write(*,*) 'i  = ', i
      write(*,*) 'j1 = ', j1
      write(*,*) 'j2 = ', j2
      write(*,*) 'nda= ', nda
      write(6, '(a)')  &
        ' *** error in check: pointer ia erroneous ***'
      Amg % ierr = AMG_ERR_IA_POINTER_INVALID
      return
    endif
    if (ja(j1).ne.i) then
      write(*,*) 'error, (ja(j1).ne.i) '
      write(*,*) 'i      = ', i
      write(*,*) 'j1     = ', j1
      write(*,*) 'ja(j1) = ', ja(j1)
      write(6, '(a)')  &
        ' *** error in check: diagonal is not stored first ***'
      Amg % ierr = AMG_ERR_DIAG_NOT_FIRST
      return
    endif
    do j = j1, j2
      i1 = ja(j)
      if (i1.lt.1 .or. i1.gt.nnu .or. icg(i1).eq.1) then
        write(*,*) 'error, (i1.lt.1 ; i1.gt.nnu ; icg(i1).eq.1) '
        write(*,*) 'i       = ', i
        write(*,*) 'nnu     = ', nnu
        write(*,*) 'i1      = ', i1
        write(*,*) 'icg(i1) = ', icg(i1)
        write(6, '(a)')  &
          ' *** error in check: pointer ja erroneous ***'
        Amg % ierr = AMG_ERR_JA_POINTER_INVALID
        return
      endif
      icg(i1) = 1
    end do
    do j = j1, j2
      icg(ja(j)) = 0
    end do
  end do
  nna = ia(nnu+1)-1

  !-----------------------------------------------------!
  !   Check for properties of a. in particular, count   !
  !   missing storage places ("new"). return if new=0   !
  !-----------------------------------------------------!
  anormm = 0.0
  anormp = 0.0
  naoff  = 0
  napos  = 0
  naneg  = 0
  nazer  = 0
  new    = 0
  do i = 1, nnu
    icg(i) = ia(i+1)-ia(i)
  end do

  do i = 1, nnu
    d = a(ia(i))
    if(d .le. 0.0) then
      write(6, '(a)')  &
        ' *** error in check: diagonal is non-positive ***'
      Amg % ierr = AMG_ERR_DIAG_NOT_POSITIVE
      return
    end if
    deps = d * 1.0e-12
    rowsum = d
    anormp = anormp+2.0*d**2
    do j = ia(i)+1, ia(i+1)-1
      rowsum = rowsum+a(j)
      if (a(j).ge.deps) naoff = naoff+1
      i1 = ja(j)
      if(i1 .le. 0) then
        ja(j) = -i1
        cycle
      end if
      if(i1 .ge. i) then
        found = .false.
        do j1 = ia(i1)+1, ia(i1+1) - 1
          if(ja(j1) .eq. i) then
            ja(j1) = -ja(j1)
            anormm = anormm+(a(j)-a(j1))**2
            anormp = anormp+(a(j)+a(j1))**2
            found = .true.
          end if
        end do
        if(found) cycle
      end if
      ja(j)=-ja(j)
      anormm = anormm+a(j)**2
      anormp = anormp+a(j)**2
      new = new+1
      icg(i1) = icg(i1)+1
    end do
    if (rowsum.gt.deps) then
      napos = napos+1
    else if (rowsum.lt.-deps) then
      naneg = naneg+1
    else
      nazer = nazer+1
    endif
  end do
  anormm = sqrt(anormm)
  anormp = sqrt(anormp)
  asym = anormm/anormp
  if (asym.le.1.0d-12) asym = 0.0

  !-------------------!
  !   Messages on a   !
  !-------------------!
# ifdef VERBOSE
    if(asym .eq. 0.0) write(6, '(a)')  &
      ' check: a probably symmetric'
    if(asym .ne. 0.0) write(6, '(a,d11.3)')  &
      ' check: a probably not symmetric. measure:', asym
    if(naoff .gt. 0)  write(6, '(a,i6,a)')  &
      ' check: a probably not pos. type:', naoff,  &
      ' off-diagonal elements positive'
    if(naneg .gt. 0)    write(6, '(a,i6,a)')  &
      ' check: a probably not pos. type:', naneg,  &
      ' rowsums negative'
    if(nazer .eq. nnu)  write(6, '(a)')  &
      ' check: a probably singular - rowsums are zero'
    if(naoff .eq. 0 .and. naneg.eq.0) write(6, '(a)')  &
      ' check: a probably positive type'
# endif

  !--------------!
  !   Warnings   !
  !--------------!
  if(Amg % isym .eq. AMG_SYMMETRIC_MATRIX .and. asym .ne. 0.0) then
    write(6, '(a)')  &
      ' --- warning 1 in check: param matrix may be bad ---'
    Amg % ierr = AMG_ERR_MATRIX_INVALID
  endif
  if(Amg % isym .eq. AMG_NON_SYMMETRIC_MATRIX .and. asym .eq. 0.0) then
    write(6, '(a)')  &
      ' --- warning 2 in check: param matrix may be bad ---'
    Amg % ierr = AMG_ERR_MATRIX_INVALID
  endif
  if(Amg % irow0 .eq. AMG_SINGULAR_MATRIX     .and. nazer .ne. nnu) then
    write(6, '(a)')  &
      ' --- warning 3 in check: param matrix may be bad ---'
    Amg % ierr = AMG_ERR_MATRIX_INVALID
  endif
  if(Amg % irow0 .eq. AMG_NON_SINGULAR_MATRIX .and. nazer .eq. nnu) then
    write(6, '(a)')  &
      ' --- warning 4 in check: param matrix may be bad ---'
    Amg % ierr = AMG_ERR_MATRIX_INVALID
  endif

  if(new .le. 0) then
#   ifdef VERBOSE
      write(6, '(a)')  &
        ' check: matrix a was symmetrically stored'
#   endif
    ! Remove pairs of zeroes
    call Amg % Truncate_Operator(1, 0, a, ia, ja)
    return
  end if

  !--------------------------------------!
  !   Replace a by symmetrized version   !
  !--------------------------------------!
  if(nna+new .ge. nda) then
    write(6, '(a)')  &
      ' *** error in check: nda too small ***'
    Amg % ierr = AMG_ERR_DIM_A_TOO_SMALL
    return
  end if
  if(nna+new .ge. nda) then
    write(6, '(a)')  &
      ' *** error in check: nda too small ***'
    Amg % ierr = AMG_ERR_DIM_JA_TOO_SMALL
    return
  end if

  !-----------------------------!
  !   Extend matrix a in situ   !
  !-----------------------------!
  ifg(1) = 1
  do i = 2, nnu + 1
    ifg(i) = ifg(i-1)+icg(i-1)
  end do
  do i = nnu, 1, -1
    ishift = ifg(i)-ia(i)
    do j = ia(i+1) - 1, ia(i), -1
      a(j+ishift)  =  a(j)
      ja(j+ishift) = ja(j)
    end do
    icg(i)  = ia(i+1)+ishift
    ia(i+1) = ifg(i+1)
  end do

  !------------------------------------------------------------------!
  !   Symmetrize matrix a:  isym=1: copy missing transpose entries   !
  !                         isym=2: fill in zeroes                   !
  !------------------------------------------------------------------!
  if(Amg % isym .eq. AMG_NON_SYMMETRIC_MATRIX) then
    do i = 1, nnu
      j2 = icg(i)-1
      do j = ia(i) + 1, j2
        i1 = ja(j)
        if(i1 .le. 0) then
          ja(j) = -i1
          jnew=icg(ja(j))
          a(jnew)= 0.0
          ja(jnew) = i
          icg(ja(j))=jnew+1
        end if
      end do
    end do
    write(6, '(a,a,i6,a)')                        &
      ' check: storage of a has been symmetrized by introducing',  &
      new,' zeroes'
  else
    write(6, '(a,i5,a,a)')                   &
      ' *** warng in check:', new, ' a-entries missing ***',  &
      ' missing transpose connections will be filled in'
    Amg % ierr = AMG_ERR_MISSING_A_ENTRY
    do i = 1, nnu
      j2 = icg(i)-1
      do j = ia(i) + 1, j2
        i1 = ja(j)
        if(i1 .le. 0) then
          ja(j) = -i1
          jnew=icg(ja(j))
          a(jnew)= 0.0
          ja(jnew) = i
          icg(ja(j))=jnew+1
        end if
      end do
    end do
  endif

  end subroutine
