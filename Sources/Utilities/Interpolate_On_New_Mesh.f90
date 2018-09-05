!======================================================================*
      program sub_ini
!----------------------------------------------------------------------*
!  This program uses solution from a previous computation to generate
!  input files that initiate velocity, pressure and temperature fields
!  for new computations. It reads files *.xyz, *.U__, *.V__, *.W__, 
!  *.P__ and *.T_ and creates *.ini files that are read by
!  LoaIni.f90.   
!----------------------------------------------------------------------*
      implicit none
!======================================================================*

  real,allocatable :: xc(:),yc(:),zc(:)
  real,allocatable :: Sx(:),Sy(:),Sz(:)
  real,allocatable :: volume(:)            ! cell's volume
  real,allocatable :: delta(:)             ! delta (max(dx,dy,dz))
  real,allocatable :: Dx(:),Dy(:),Dz(:)
  real,allocatable :: xsp(:),ysp(:),zsp(:) ! face coordinates
  real,allocatable :: WallDs(:), f(:)
  real             :: Xmax, Xmin, Ymin, Ymax, Zmin, Zmax

  integer   :: NC, NS, ND, NN                    ! num. of nodes and cells
  integer   :: NbC, NSsh, Nmat

  integer,allocatable :: material(:)     ! material markers
  integer,allocatable :: SideC(:,:)      !  c0, c1, c2

  integer          :: l1, n, IND
  integer          :: j, k, c, c1, c2, s
  integer          :: NCold
  character*80     :: nameIn
  character*80     :: namSav
  real,allocatable :: Xold(:),Yold(:),Zold(:)
  real,allocatable :: Uold(:),Vold(:),Wold(:),Told(:)
  real,allocatable :: UCold(:),VCold(:),WCold(:),TCold(:)
  real,allocatable :: UCoold(:),VCoold(:),WCoold(:),TCoold(:)
  real,allocatable :: Uoold(:),Voold(:),Woold(:),Toold(:)
  real,allocatable :: UDoold(:),VDoold(:),WDoold(:),TDoold(:)
  real,allocatable :: UXold(:),VXold(:),WXold(:),TXold(:)
  real,allocatable :: UXoold(:),VXoold(:),WXoold(:),TXoold(:)
  real,allocatable :: Pold(:) 
  real,allocatable :: PPold(:)
  real,allocatable :: Pxold(:),Pyold(:),Pzold(:)
!---- Variables for ReadC:
  character  :: namU*80, namV*80, namW*80, naOut*80, naIn*80, namFin*80
  character  :: answer*80, nameOut*80, namP*80, namT*80, name*80
!----------------------------------------------------------------------!
! The answer name is case dependent
!----------------------------------------------------------------------!

  print *,' Enter name of case which solution you interpolate on new mesh (without ext.) : '
  read(*,*) name 
  print *,' Enter name of present case (without ext.): '
  read(*,*) namSav 
  print *,' Enter number of subdomains: '
  read(*,*) ND 
  print *,' Enter 1 if temperature field is to be interpolated and 0 if not : '
  read(*,*) IND

  answer = name
  naOut = answer
  call Name_File(name, naOut, '.xyz')

  answer = name
  namU = answer
  namU(len_trim(answer)+1:len_trim(answer)+4)='.U__'

  answer = name
  namV = answer
  namV(len_trim(answer)+1:len_trim(answer)+4)='.V__'

  answer = name
  namW = answer
  namW(len_trim(answer)+1:len_trim(answer)+4)='.W__'

  answer = name
  namP = answer
  namP(len_trim(answer)+1:len_trim(answer)+4)='.P__'

  answer = name
  namT = answer
  namT(len_trim(answer)+1:len_trim(answer)+4)='.T__'

  open(5, file=naOut)
  read(5,*) NCold

  allocate (Xold(NCold)); Xold = 0.0
  allocate (Yold(NCold)); Yold = 0.0
  allocate (Zold(NCold)); Zold = 0.0
  allocate (Uold(NCold)); Uold = 0.0
  allocate (Vold(NCold)); Vold = 0.0
  allocate (Wold(NCold)); Wold = 0.0
  allocate (Told(NCold)); Told = 0.0
  allocate (Uoold(NCold)); Uoold = 0.0
  allocate (Voold(NCold)); Voold = 0.0
  allocate (Woold(NCold)); Woold = 0.0
  allocate (Toold(NCold)); Toold = 0.0
  allocate (UDoold(NCold)); UDoold = 0.0
  allocate (VDoold(NCold)); VDoold = 0.0
  allocate (WDoold(NCold)); WDoold = 0.0
  allocate (TDoold(NCold)); TDoold = 0.0
  allocate (UCold(NCold)); UCold = 0.0
  allocate (VCold(NCold)); VCold = 0.0
  allocate (WCold(NCold)); WCold = 0.0
  allocate (TCold(NCold)); TCold = 0.0
  allocate (UCoold(NCold)); UCoold = 0.0
  allocate (VCoold(NCold)); VCoold = 0.0
  allocate (WCoold(NCold)); WCoold = 0.0
  allocate (TCoold(NCold)); TCoold = 0.0
  allocate (UXold(NCold)); UXold = 0.0
  allocate (VXold(NCold)); VXold = 0.0
  allocate (WXold(NCold)); WXold = 0.0
  allocate (TXold(NCold)); TXold = 0.0
  allocate (UXoold(NCold)); UXoold = 0.0
  allocate (VXoold(NCold)); VXoold = 0.0
  allocate (WXoold(NCold)); WXoold = 0.0
  allocate (TXoold(NCold)); TXoold = 0.0
  allocate (Pold(NCold)); Pold = 0.0
  allocate (PPold(NCold)); PPold = 0.0
  allocate (Pxold(NCold)); Pxold = 0.0
  allocate (Pyold(NCold)); Pyold = 0.0
  allocate (Pzold(NCold)); Pzold = 0.0
  n = 0
  j = NCold
  do k = 1, j
    if(mod(k,500000) .eq. 0) print *, (100.*k/(1.*j)), '% complete...'
    read(5,*) Xold(k), Yold(k), Zold(k)
    n = n + 1
  end do
  close(5)
  print *, ' Finished with reading case.xyz file', n
  open(5, file=namU)
  do k = 1, j
    if(mod(k,500000) .eq. 0) print *, (100.*k/(1.*j)), '% complete...'
    read(5,*) Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k)
  end do
  close(5)
  print *, ' Finished with reading case.U__ file'
  open(5, file=namV)
  do k = 1, j
    if(mod(k,500000) .eq. 0) print *, (100.*k/(1.*j)), '% complete...'
    read(5,*) Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k)
  end do
  close(5)
  print *, ' Finished with reading case.V__ file'
  open(5, file=namW)
  do k = 1, j
    if(mod(k,500000) .eq. 0) print *, (100.*k/(1.*j)), '% complete...'
      read(5,*) Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k)
    end do
  close(5)
  print *, ' Finished with reading case.W__ file'
  open(5, file=namP)
  do k = 1, j
    if(mod(k,500000) .eq. 0) print *, (100.*k/(1.*j)), '% complete...'
      read(5,*) Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
    end do
  close(5)
  print *, ' Finished with reading case.P__ file'

  if(IND .eq. 1) then
    open(5, file=namT)
    do k = 1, j
      if(mod(k,500000) .eq. 0) print *, (100.*k/(1.*j)), '% complete...'
        read(5,*) Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k) 
      end do
    close(5)
    print *, ' Finished with reading case.T__ file'
    print *, 'LoaInI: finished with reading the files'
  end if

  NameOut = namSav
  if(ND > 1) then
    NameOut(len_trim(namSav)+1:len_trim(namSav)+5)="-0000"
    l1=len_trim(NameOut)
  end if 

  do j = 1, ND 
    if(ND > 1) then
      if(j  <  10) then
          write(NameOut(l1  :l1),'(I1)') j
      else if(j  < 100) then
        write(NameOut(l1-1:l1),'(I2)') j
      else
        write(NameOut(l1-2:l1),'(I3)') j
      end if
    end if  
    print *, NameOut

    name = NameOut
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
!     Read the binary file with the     *
!       connections between cells       *
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
    nameIn=name
    nameIn(len_trim(name)+1:len_trim(name)+4)='.cns'
    open(9, file=nameIn,FORM='unformatted')
    print *, '# Now reading the binary .cns file: ', nameIn

!///// number of cells, boundary cells and sides
    read(9) NC
    read(9) NbC
    read(9) NS
    read(9) NSsh
    read(9) Nmat

!///// cell materials
    allocate (material(-NbC:NC))
    allocate (SideC(0:2,NS))
    read(9) (material(c), c=1,NC)
    read(9) (material(c), c=-1,-NBC,-1)

    close(9)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
!     Read the binary file with     *
!       geometrical quantities      *
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
    nameIn = name
    nameIn(len_trim(name)+1:len_trim(name)+4)='.geo'
    open(9, file=nameIn, FORM='unformatted')
    print *, '# Now reading the binary .geo file: ', nameIn

    allocate (xc(-NbC:NC))
    allocate (yc(-NbC:NC))
    allocate (zc(-NbC:NC))
    allocate (volume(-NbC:NC))
    allocate (delta(-NbC:NC))
    allocate (WallDs(NS))
    allocate (Sx(NS))
    allocate (Sy(NS))
    allocate (Sz(NS))
    allocate (Dx(NS))
    allocate (Dy(NS))
    allocate (Dz(NS))
    allocate (f(NS))
    allocate (xsp(NS))
    allocate (ysp(NS))
    allocate (zsp(NS))

    read(9) (xc(c), c=1,NC)
    read(9) (yc(c), c=1,NC)
    read(9) (zc(c), c=1,NC)
    read(9) (xc(c), c=-1,-NBC,-1)
    read(9) (yc(c), c=-1,-NBC,-1)
    read(9) (zc(c), c=-1,-NBC,-1)

    read(9) (volume(c), c=1,NC)
    read(9) (delta(c),  c=1,NC)

    read(9) (WallDs(c), c=1,NC)

    read(9) (Sx(s), s=1,NS)
    read(9) (Sy(s), s=1,NS)
    read(9) (Sz(s), s=1,NS)

    read(9) (Dx(s), s=1,NS)
    read(9) (Dy(s), s=1,NS)
    read(9) (Dz(s), s=1,NS)

    read(9) (f(s), s=1,NS)

    read(9) (xsp(s), s=1,NS)
    read(9) (ysp(s), s=1,NS)
    read(9) (zsp(s), s=1,NS)

    close(9)
 
    Xmax = -1.0e+6
    Xmin = 1.0e+6
    Ymax = -1.0e+6
    Ymin = 1.0e+6
    Zmax = -1.0e+6
    Zmin = 1.0e+6

    do s = 1, NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      Xmax = max(xsp(s),Xmax)
      Xmin = min(xsp(s),Xmin)
      Ymax = max(ysp(s),Ymax)
      Ymin = min(ysp(s),Ymin)
      Zmax = max(zsp(s),Zmax)
      Zmin = min(zsp(s),Zmin)
    end do  
    
    namFin=name
    namFin(len_trim(name)+1:len_trim(name)+4)='.ini'
    print *, namFin
    open(9,file = namFin)
    n = 0
    k = 0
    do k = 1, NCold
      if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
        if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
          if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
            n = n + 1
          end if
        end if
      end if
    end do

    write(9,*) n
    NN = 0
    if(IND .eq. 1) then
      do k = 1, NCold
        if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
          if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
            if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
              NN = NN + 1
              write(9,*) Xold(k), Yold(k), Zold(k), &
                                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                                  Told(k), Toold(k), TCold(k), TCoold(k), TDoold(k), TXold(k), TXoold(k), &
                                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
            end if
          end if
        end if
      end do
    else if(IND .eq. 0) then
      do k = 1, NCold
        if(Xold(k) <= Xmax.and.Xold(k) >= Xmin) then
          if(Yold(k) <= Ymax.and.Yold(k) >= Ymin) then
            if(Zold(k) <= Zmax.and.Zold(k) >= Zmin) then
              NN = NN + 1
              write(9,*) Xold(k), Yold(k), Zold(k), &
                                  Uold(k), Uoold(k), UCold(k), UCoold(k), UDoold(k), UXold(k), UXoold(k), &
                                  Vold(k), Voold(k), VCold(k), VCoold(k), VDoold(k), VXold(k), VXoold(k), &
                                  Wold(k), Woold(k), WCold(k), WCoold(k), WDoold(k), WXold(k), WXoold(k), &
                                  Pold(k), PPold(k), Pxold(k), Pyold(k), Pzold(k)
            end if
          end if
        end if
      end do
    end if
    print *, 'The number of cells in subdomain is: ', NC, 'Found cells: ', NN
    close(9)
    deallocate(material)
    deallocate (xc)
    deallocate (yc)
    deallocate (zc)
    deallocate (volume)
    deallocate (delta)
    deallocate (WallDs)
    deallocate (Sx)
    deallocate (Sy)
    deallocate (Sz)
    deallocate (Dx)
    deallocate (Dy)
    deallocate (Dz)
    deallocate (f)
    deallocate (xsp)
    deallocate (ysp)
    deallocate (zsp)
    deallocate (SideC)
  end do

  deallocate(Xold)
  deallocate(Yold)
  deallocate(Zold)
  deallocate(Uold)
  deallocate(Vold)
  deallocate(Wold)
  deallocate(Told)
  deallocate(Uoold)
  deallocate(Voold)
  deallocate(Woold)
  deallocate(Toold)
  deallocate(UDoold)
  deallocate(VDoold)
  deallocate(WDoold)
  deallocate(TDoold)
  deallocate(UCold)
  deallocate(VCold)
  deallocate(WCold)
  deallocate(TCold)
  deallocate(UCoold)
  deallocate(VCoold)
  deallocate(WCoold)
  deallocate(TCoold)
  deallocate(UXold)
  deallocate(VXold)
  deallocate(WXold)
  deallocate(TXold)
  deallocate(UXoold)
  deallocate(VXoold)
  deallocate(WXoold)
  deallocate(TXoold)
  deallocate(Pold)
  deallocate(PPold)
  deallocate(Pxold)
  deallocate(Pyold)
  deallocate(Pzold)

  end program

!======================================================================!
  subroutine Name_File(name, namOut, ext, lext)
!----------------------------------------------------------------------!
!   Creates the file name depending on the subdomain and file type.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
!  use all_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer       :: sub, lext
  character*(*) :: ext
  character*(*) :: namOut
  character     :: name*80
!-------------------------------[Locals]-------------------------------!
  integer   :: c
  character :: numb*4
!======================================================================!

  namOut = name

!  if(sub .eq. 0) then
    namOut(len_trim(name)+1:len_trim(name)+1+lext-1) = ext(1:lext)
!  else
!    write(numb,'(I4)') sub
!    write(namOut(len_trim(name)+1:len_trim(name)+5),'(A5)') '-0000'
!    do c=1,4
!      if( numb(c:c) >= '0' .and. numb(c:c) <= '9' )                 &
!        namOut(len_trim(name)+1+c:len_trim(name)+1+c) = numb(c:c)
!    end do
!    namOut(len_trim(name)+6:len_trim(name)+6+lext-1) = ext(1:lext)
!  end if

  end subroutine

