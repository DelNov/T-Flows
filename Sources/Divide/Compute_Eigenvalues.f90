!==============================================================================!
  subroutine Compute_Eigenvalues(a,n,np,d,v,nrot) 
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: NMAX=500
  integer            :: n,np,nrot
  integer            :: i,ip,iq,j 
  real               :: a(np,np),d(np),v(np,np) 
  real               :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX) 
!------------------------------------------------------------------------------!
!   Uses Jacobi algorithm to computes eigenvalues and eigenvectors of a real   !
!   symmetric matrix a, which is of size n by n, stored in a physical np by    !
!   np array. On output, elements of a above the diagonal are destroyed.       !
!   d returns the eigenvalues of a in its  rst n elements. v is a matrix with  !
!   the same logical and physical dimensions as a, whose columns contain, on   !
!   output, the normalized eigenvectors of a. nrot returns the number of       !
!   Jacobi rotations that were required.                                       !
!==============================================================================!

  do ip=1,n     ! initialize to the identity matrix. 
    do iq=1,n 
      v(ip,iq)=0. 
    end do 
    v(ip,ip)=1. 
  end do 

  do ip=1,n 
    b(ip)=a(ip,ip) ! initialize b and d to the diagonal of a. 
    d(ip)=b(ip) 
    z(ip)=0.       ! this vector will accumulate terms of the form 
                   ! tapq as in equation (11.1.14). 
  end do 

  nrot=0 
  do i=1,50 

    sm=0. 
    do ip=1,n-1  ! sum o -diagonal elements. 
      do iq=ip+1,n 
        sm=sm+abs(a(ip,iq)) 
      end do 
    end do 
    if(sm.eq.0.) return  ! the normal return, which relies on quadratic 
                         ! convergence to machine under ow. 
    if(i.lt.4) then 
      tresh=0.2*sm/n**2  ! ...on the  rst three sweeps. 
    else 
      tresh=0.           ! ...thereafter. 
    endif 

    do ip=1,n-1 
      do iq=ip+1,n 
        g=100.*abs(a(ip,iq))     ! after four sweeps, skip the rotation if the 
                                 ! offdiagonal element is small. 
        if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) &
           .and.(abs(d(iq))+g.eq.abs(d(iq)))) then 
          a(ip,iq)=0. 
        else if(abs(a(ip,iq)).gt.tresh) then 
          h=d(iq)-d(ip) 
          if(abs(h)+g.eq.abs(h)) then
            t=a(ip,iq)/h         ! t = 1/(2 theta) 
          else 
            theta=0.5*h/a(ip,iq) ! equation (11.1.10). 
            t=1./(abs(theta)+sqrt(1.+theta**2)) 
            if(theta.lt.0.) t=-t 
          endif 
          c=1./sqrt(1+t**2) 
          s=t*c 
          tau=s/(1.+c) 
          h=t*a(ip,iq) 
          z(ip)=z(ip)-h 
          z(iq)=z(iq)+h 
          d(ip)=d(ip)-h 
          d(iq)=d(iq)+h 
          a(ip,iq)=0. 
          do j=1,ip-1     ! case of rotations 1   j <p. 
            g=a(j,ip) 
            h=a(j,iq) 
            a(j,ip)=g-s*(h+g*tau) 
            a(j,iq)=h+s*(g-h*tau) 
          end do 
          do j=ip+1,iq-1  ! case of rotations p <j<q. 
            g=a(ip,j) 
            h=a(j,iq) 
            a(ip,j)=g-s*(h+g*tau) 
            a(j,iq)=h+s*(g-h*tau) 
          end do 
          do j=iq+1,n     ! case of rotations q <j   n. 
            g=a(ip,j) 
            h=a(iq,j) 
            a(ip,j)=g-s*(h+g*tau) 
            a(iq,j)=h+s*(g-h*tau) 
          end do 
          do j=1,n 
            g=v(j,ip) 
            h=v(j,iq) 
            v(j,ip)=g-s*(h+g*tau) 
            v(j,iq)=h+s*(g-h*tau) 
          end do 
          nrot=nrot+1 
        endif 
      end do 
    end do 

    do ip=1,n 
      b(ip)=b(ip)+z(ip) 
      d(ip)=b(ip)          ! update d with the sum of tapq, z(ip)=0. 
                           ! and reinitialize z. 
    end do 
  end do 

  print *, '# Too many Jacobi iterations in Compute_Eigenvalues' 

  return 

  end subroutine Compute_Eigenvalues
