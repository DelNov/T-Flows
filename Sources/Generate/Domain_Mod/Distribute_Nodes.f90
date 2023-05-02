!==============================================================================!
  subroutine Distribute_Nodes(Dom, Grid, b, w,   &
                              is, js, ks, ie, je, ke)
!------------------------------------------------------------------------------!
!   Places the nodes on the line defined with local block position             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom
  type(Grid_Type)     :: Grid
  integer, intent(in) :: b, is, js, ks, ie, je, ke
  real,    intent(in) :: w
!----------------------------------[Calling]-----------------------------------!
  real :: atanh
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, ni, nj, nk, i, j, k, node, case
  real    :: x0, y0, z0, delx, dely, delz, t, dt, ddt, pr, xi
!==============================================================================!

  ni = Dom % blocks(b) % resolutions(1)
  nj = Dom % blocks(b) % resolutions(2)
  nk = Dom % blocks(b) % resolutions(3)

  x0   = Grid % xn(Grid % n_nodes+(ks-1)*ni*nj+(js-1)*ni+is)
  y0   = Grid % yn(Grid % n_nodes+(ks-1)*ni*nj+(js-1)*ni+is)
  z0   = Grid % zn(Grid % n_nodes+(ks-1)*ni*nj+(js-1)*ni+is)
  delx = Grid % xn(Grid % n_nodes+(ke-1)*ni*nj+(je-1)*ni+ie) - x0
  dely = Grid % yn(Grid % n_nodes+(ke-1)*ni*nj+(je-1)*ni+ie) - y0
  delz = Grid % zn(Grid % n_nodes+(ke-1)*ni*nj+(je-1)*ni+ie) - z0

  n = max( (ie-is), (je-js),  (ke-ks) )

  if(n < 2) return

  !-------------------------!
  !   Linear distribution   !
  !-------------------------!
  if(w > 0.0) then
    ddt = ( 2.0*(1.0-w) ) / ( real(n)*(real(n)-1.0)*(1.0+w) )
    t=0.0
    do i=is,ie
      do j=js,je
        do k=ks,ke
          if( ie .ne. is ) then
            dt=1.0/real(n)+(real(i)-0.5*(real(n)+1)) * ddt
            t=t+dt
            node = Grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i+1
            if( (i  < ie).and.(Grid % xn(node) .gt. TERA) ) then
              Grid % xn(node) = x0 + t*delx
              Grid % yn(node) = y0 + t*dely
              Grid % zn(node) = z0 + t*delz
            end if
          end if
          if( je .ne. js ) then
            dt=1.0/real(n)+(real(j)-0.5*(real(n)+1)) * ddt
            t=t+dt
            node = Grid % n_nodes + (k-1)*ni*nj + (j-0)*ni + i
            if( (j  < je).and.(Grid % xn(node) .gt. TERA) ) then
              Grid % xn(node) = x0 + t*delx
              Grid % yn(node) = y0 + t*dely
              Grid % zn(node) = z0 + t*delz
            end if
          end if
          if( ke .ne. ks ) then
            dt=1.0/real(n)+(real(k)-0.5*(real(n)+1)) * ddt
            t=t+dt
            node = Grid % n_nodes + (k-0)*ni*nj + (j-1)*ni + i
            if( (k  < ke).and.(Grid % xn(node) .gt. TERA) ) then
              Grid % xn(node) = x0 + t*delx
              Grid % yn(node) = y0 + t*dely
              Grid % zn(node) = z0 + t*delz
            end if
          end if
        end do
      end do
    end do

  !-----------------------------!
  !   Hyperbolic distribution   !
  !-----------------------------!
  else
    case = 0
    if     ((w  >  -0.5).and.(w <=  -0.25)) then
      pr = 1.0 - abs(0.5 - abs(w))
      case = 1
    else if((w >=  -0.75).and.(w  <  -0.5)) then
      pr = 1.0 - abs(0.5 - abs(w))
      case = 2
    else
      pr = -w
      case = 3
    end if

    do i=is,ie
      do j=js,je
        do k=ks,ke
          if( ie .ne. is ) then
            if(case .eq. 1) xi =       - 1.0*real(i)/real(n)
            if(case .eq. 2) xi =   1.0 - 1.0*real(i)/real(n)
            if(case .eq. 3) xi = - 1.0 + 2.0*real(i)/real(n)
            node = Grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i+1
            if( (i  < ie).and.(Grid % xn(node) .gt. TERA) ) then
              if    (case .eq. 1) then
                Grid % xn(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 2) then
                Grid % xn(node) = x0  &
                                + delx - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + dely - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + delz - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 3) then
                Grid % xn(node) = x0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
              end if
            end if
          end if
          if( je .ne. js ) then
            if(case .eq. 1) xi =        - 1.0*real(j)/real(n)
            if(case .eq. 2) xi =    1.0 - 1.0*real(j)/real(n)
            if(case .eq. 3) xi = -  1.0 + 2.0*real(j)/real(n)
            node = Grid % n_nodes + (k-1)*ni*nj + (j-0)*ni + i
            if( (j  < je).and.(Grid % xn(node) .gt. TERA) ) then
              if    (case .eq. 1) then
                Grid % xn(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 2) then
                Grid % xn(node) = x0  &
                                + delx - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + dely - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + delz - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 3) then
                Grid % xn(node) = x0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
              end if
            end if
          end if
          if( ke .ne. ks ) then
            if(case .eq. 1) xi =       - 1.0*real(k)/real(n)
            if(case .eq. 2) xi =   1.0 - 1.0*real(k)/real(n)
            if(case .eq. 3) xi = - 1.0 + 2.0*real(k)/real(n)
            node = Grid % n_nodes + (k-0)*ni*nj + (j-1)*ni + i
            if( (k  < ke).and.(Grid % xn(node) .gt. TERA) ) then
              if    (case .eq. 1) then
                Grid % xn(node) = x0 - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0 - (tanh(xi*atanh(pr))/pr)*dely
                Grid % yn(node) = z0 - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 2) then
                Grid % xn(node) = x0  &
                                + delx - (tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + dely - (tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + delz - (tanh(xi*atanh(pr))/pr)*delz
              elseif(case .eq. 3) then
                Grid % xn(node) = x0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
                Grid % yn(node) = y0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*dely
                Grid % zn(node) = z0  &
                                + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delz
              end if
            end if
          end if
        end do
      end do
    end do

  end if

  end subroutine
