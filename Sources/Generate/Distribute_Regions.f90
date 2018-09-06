!==============================================================================!
  subroutine Distribute_Regions(dom, grid)
!------------------------------------------------------------------------------!
!   Distribute regions (defined in .dom file) to boundary conditions           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Gen_Mod
  use Domain_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, i, j, k, n, c, r
  integer :: n_mat, n_bnd         ! number of materials and boundary conditions
  integer :: is, js, ks, ie, je, ke, face 
  integer :: ci, cj, ck
  logical :: found
!==============================================================================!

  !-----------------------------------------!
  !   Insertion of the boundary condition   ! 
  !        and materials information        !
  !-----------------------------------------!

  ! This is too much memory but that's OK 
  !  (+1 is to store the default values)
  allocate(grid % bnd_cond % name(dom % n_regions + 1))

  ! Set the bare bones - one material minimum, and call it "AIR"
  n_bnd = 0
  n_mat = 1
  grid % material % name = "AIR"

  do n = 1, dom % n_regions

    b = dom % regions(n) % block

    ! Block resolution
    ci = dom % blocks(b) % resolutions(1)-1
    cj = dom % blocks(b) % resolutions(2)-1
    ck = dom % blocks(b) % resolutions(3)-1

    ! Default values
    is = 1
    ie = ci
    js = 1
    je = cj
    ks = 1
    ke = ck

    ! Boundary conditions prescribed with mnemonics
    if(dom % regions(n) % face .eq. 'IMIN') then
      ie=1 
      face = 5
    else if(dom % regions(n) % face .eq. 'IMAX') then 
      is=ci
      face = 3
    else if(dom % regions(n) % face .eq. 'JMIN') then 
      je=1
      face = 2
    else if(dom % regions(n) % face .eq. 'JMAX') then 
      js=cj
      face = 4
    else if(dom % regions(n) % face .eq. 'KMIN') then 
      ke=1
      face = 1
    else if(dom % regions(n) % face .eq. 'KMAX') then 
      ks=ck
      face = 6

    ! Boundary conditions  prescribed explicitly
    !  (error prone and difficult, but might be usefull)
    else   
      is = dom % regions(n) % is
      js = dom % regions(n) % js
      ks = dom % regions(n) % ks
      ie = dom % regions(n) % ie
      je = dom % regions(n) % je
      ke = dom % regions(n) % ke
      face = 0
      if( (is .eq. ie).and.(is .eq.  1) ) face=5
      if( (is .eq. ie).and.(is .eq. ci) ) face=3
      if( (js .eq. je).and.(js .eq.  1) ) face=2
      if( (js .eq. je).and.(js .eq. cj) ) face=4
      if( (ks .eq. ke).and.(ks .eq.  1) ) face=1
      if( (ks .eq. ke).and.(ks .eq. ck) ) face=6
    end if

    ! Store boundary condition 
    if(face .ne. 0) then  

      found = .false. 
      do r=1,n_bnd
        if( grid % bnd_cond % name(r) .eq.   &
            dom % regions(n) % name ) found = .true.
      end do
      if( .not. found) then
        n_bnd = n_bnd + 1
        grid % bnd_cond % name(n_bnd) = dom % regions(n) % name
      end if

      do i=is,ie
        do j=js,je
          do k=ks,ke
            c = dom % blocks(b) % n_cells + (k-1)*ci*cj + (j-1)*ci + i   
            grid % cells_c(face,c) = -n_bnd
          end do
        end do
      end do

    end if

  end do  !  n_regions

  ! Store the number of boundary conditions
  grid % n_bnd_cond  = n_bnd

  print *, '#==================================================='
  print *, '# Found following boundary conditions:'
  print *, '#---------------------------------------------------'
  do n = 1, grid % n_bnd_cond
    print *, '# ', trim(grid % bnd_cond % name(n))
  end do
  print *, '#---------------------------------------------------'

  end subroutine
