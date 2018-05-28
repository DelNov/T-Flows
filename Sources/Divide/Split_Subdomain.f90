!==============================================================================!
  subroutine Split_Subdomain(grid, sub, n_parts)
!------------------------------------------------------------------------------!
!   Splits a (sub) domain by a selected technique.                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use gen_mod 
  use Div_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub      ! subdomain to be splitted
  integer         :: n_parts  ! number of new partitions
!-----------------------------------[Locals]-----------------------------------!
  character :: dir                         ! direction for splitting
  integer   :: c, ic, j
  integer   :: n_filled
  real      :: xmax,ymax,zmax,xmin,ymin,zmin, delx,dely,delz,dxyz 
!==============================================================================!

  !-----------------------------------------! 
  !   Find the smallest moment of inertia   !
  !-----------------------------------------! 
  if(division_algorithm .eq. INERTIAL) then
    call Inertia(grid, sub)
  end if

  !----------------------------------------------! 
  !   Find the largest dimension of the domain   !
  !----------------------------------------------! 
  if(division_algorithm .eq. COORDINATE) then
    xmax=-HUGE
    ymax=-HUGE
    zmax=-HUGE
    xmin=+HUGE
    ymin=+HUGE
    zmin=+HUGE
    do c = 1, grid % n_cells
      if(proces(c) .eq. sub) then
        xmax=max(xmax, grid % xc(c))
        ymax=max(ymax, grid % yc(c))
        zmax=max(zmax, grid % zc(c))
        xmin=min(xmin, grid % xc(c))
        ymin=min(ymin, grid % yc(c))
        zmin=min(zmin, grid % zc(c))
      end if
    end do
    delx = xmax - xmin
    dely = ymax - ymin
    delz = zmax - zmin

    dxyz = max(delx,dely,delz)
    if(delz .eq. dxyz) dir='z'
    if(dely .eq. dxyz) dir='y'
    if(delx .eq. dxyz) dir='x'
  end if

  do j = 1, n_parts-1

    n_filled   = 0

    !------------------------!
    !   Fill the subdomain   !
    !------------------------!
    if(division_algorithm .eq. COORDINATE) then

      if(dir .eq. 'x') then
        do c = 1, grid % n_cells
          ic=ix(c)
          if(proces(ic) .eq. sub) then
            proces(ic) = n_sub+j
            n_filled = n_filled + 1
            if(n_filled >= sub_n_cells(sub)/(n_parts-j+1)) goto 2 
          end if
        end do
      end if

      if(dir .eq. 'y') then
        do c = 1, grid % n_cells
          ic=iy(c)
          if(proces(ic) .eq. sub) then
            proces(ic) = n_sub+j
            n_filled = n_filled + 1
            if(n_filled >= sub_n_cells(sub)/(n_parts-j+1)) goto 2 
          end if
        end do
      end if

      if(dir .eq. 'z') then
        do c = 1, grid % n_cells
          ic=iz(c)
          if(proces(ic) .eq. sub) then
            proces(ic) = n_sub+j
            n_filled = n_filled + 1
            if(n_filled >= sub_n_cells(sub)/(n_parts-j+1)) goto 2 
          end if
        end do
      end if
    end if

    if(division_algorithm .eq. INERTIAL) then
      do c = 1, grid % n_cells
        ic=iin(c)
        if(proces(ic) .eq. sub) then
          proces(ic) = n_sub+j
          n_filled = n_filled + 1
          if(n_filled >= sub_n_cells(sub)/(n_parts-j+1)) goto 2 
        end if
      end do
    end if

    !----------------------------! 
    !   Subdomain is filled up   !
    !----------------------------! 
 2  sub_n_cells(n_sub+j) = n_filled
    sub_n_cells(sub) = sub_n_cells(sub) - n_filled

  end do  ! j

  n_sub = n_sub + n_parts - 1

  end subroutine
