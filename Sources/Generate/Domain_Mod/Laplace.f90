!==============================================================================!
  subroutine Domain_Mod_Laplace(dom, grid, b,i,j,k, wx16, wx24, wx35,  &
                                                    wy16, wy24, wy35,  &
                                                    wz16, wz24, wz35)
!------------------------------------------------------------------------------!
!   Places the nodes inside the block using Laplace-like function              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: grid
  integer           :: b, i, j, k
  real              :: wx16, wx24, wx35, wy16, wy24, wy35, wz16, wz24, wz35 
!-----------------------------------[Locals]-----------------------------------!
  real    :: xt(8), yt(8), zt(8)
  integer :: ni, nj, nk, n, n1, n2, n3, n4, n5, n6
  real    :: xf1, yf1, zf1, xf2, yf2, zf2, xf3, yf3, zf3 
  real    :: xf4, yf4, zf4, xf5, yf5, zf5, xf6, yf6, zf6
!==============================================================================!

  ni = dom % blocks(b) % resolutions(1)
  nj = dom % blocks(b) % resolutions(2)
  nk = dom % blocks(b) % resolutions(3)

  do n=1,8
    xt(n) = dom % points( dom % blocks(b) % corners(n) ) % x
    yt(n) = dom % points( dom % blocks(b) % corners(n) ) % y
    zt(n) = dom % points( dom % blocks(b) % corners(n) ) % z
  end do  

  n = grid % n_nodes + (k-1)*ni*nj + (j-1)*ni + i

  ! Node numbers at the block faces
  n1 = grid % n_nodes + ( 1-1)*ni*nj + ( j-1)*ni + i     !  ->  k .eq. 1
  n2 = grid % n_nodes + ( k-1)*ni*nj + ( j-1)*ni + 1     !  ->  i .eq. 1
  n3 = grid % n_nodes + ( k-1)*ni*nj + ( 1-1)*ni + i     !  ->  j .eq. 1
  n4 = grid % n_nodes + ( k-1)*ni*nj + ( j-1)*ni + ni    !  ->  i .eq. ni
  n5 = grid % n_nodes + ( k-1)*ni*nj + (nj-1)*ni + i     !  ->  j .eq. nj
  n6 = grid % n_nodes + (nk-1)*ni*nj + ( j-1)*ni + i     !  ->  k .eq. nk

  ! Face I
  if(grid % xn(n1) .gt. TERA) then
    xf1=( (real(ni-i)*xt(1) + real(i-1)*xt(2)) * real(nj-j) +   &
          (real(ni-i)*xt(3) + real(i-1)*xt(4)) * real(j-1)  )   &
       /  (real(ni-1)*real(nj-1))
    yf1=( (real(ni-i)*yt(1) + real(i-1)*yt(2)) * real(nj-j) +   &
          (real(ni-i)*yt(3) + real(i-1)*yt(4)) * real(j-1)  )   &
       /  (real(ni-1)*real(nj-1))
    zf1=( (real(ni-i)*zt(1) + real(i-1)*zt(2)) * real(nj-j) +   &
          (real(ni-i)*zt(3) + real(i-1)*zt(4)) * real(j-1)  )   &
       /  (real(ni-1)*real(nj-1))

  else
    xf1 = grid % xn(n1)
    yf1 = grid % yn(n1)
    zf1 = grid % zn(n1)
  end if

  ! Face VI
  if(grid % xn(n6) .gt. TERA) then
    xf6=( (real(ni-i)*xt(5) + real(i-1)*xt(6)) * real(nj-j) +   &
          (real(ni-i)*xt(7) + real(i-1)*xt(8)) * real(j-1)  )   &
       /  (real(ni-1)*real(nj-1))
    yf6=( (real(ni-i)*yt(5) + real(i-1)*yt(6)) * real(nj-j) +   &
          (real(ni-i)*yt(7) + real(i-1)*yt(8)) * real(j-1)  )   &
       /  (real(ni-1)*real(nj-1))
    zf6=( (real(ni-i)*zt(5) + real(i-1)*zt(6)) * real(nj-j) +   &
          (real(ni-i)*zt(7) + real(i-1)*zt(8)) * real(j-1)  )   &
       /  (real(ni-1)*real(nj-1))
  else
    xf6 = grid % xn(n6)
    yf6 = grid % yn(n6)
    zf6 = grid % zn(n6)
  end if

  ! Face III
  if(grid % xn(n3) .gt. TERA) then
    xf3=( (real(ni-i)*xt(1) + real(i-1)*xt(2)) * real(nk-k) +   &
          (real(ni-i)*xt(5) + real(i-1)*xt(6)) * real(k-1)  )   &
       /  (real(ni-1)*real(nk-1))
    yf3=( (real(ni-i)*yt(1) + real(i-1)*yt(2)) * real(nk-k) +   &
          (real(ni-i)*yt(5) + real(i-1)*yt(6)) * real(k-1)  )   &
       /  (real(ni-1)*real(nk-1))
    zf3=( (real(ni-i)*zt(1) + real(i-1)*zt(2)) * real(nk-k) +   &
          (real(ni-i)*zt(5) + real(i-1)*zt(6)) * real(k-1)  )   &
       /  (real(ni-1)*real(nk-1))
  else
    xf3 = grid % xn(n3)
    yf3 = grid % yn(n3)
    zf3 = grid % zn(n3)
  end if

  ! Face V
  if(grid % xn(n5) .gt. TERA) then
    xf5=( (real(ni-i)*xt(3) + real(i-1)*xt(4)) * real(nk-k) +   &
          (real(ni-i)*xt(7) + real(i-1)*xt(8)) * real(k-1)  )   &
       /  (real(ni-1)*real(nk-1))
    yf5=( (real(ni-i)*yt(3) + real(i-1)*yt(4)) * real(nk-k) +   &
          (real(ni-i)*yt(7) + real(i-1)*yt(8)) * real(k-1)  )   &
       /  (real(ni-1)*real(nk-1))
    zf5=( (real(ni-i)*zt(3) + real(i-1)*zt(4)) * real(nk-k) +   &
          (real(ni-i)*zt(7) + real(i-1)*zt(8)) * real(k-1)  )   &
       /  (real(ni-1)*real(nk-1))
  else
    xf5 = grid % xn(n5)
    yf5 = grid % yn(n5)
    zf5 = grid % zn(n5)
  end if

  ! Face II
  if(grid % xn(n2) .gt. TERA) then
    xf2=( (real(nj-j)*xt(1) + real(j-1)*xt(3)) * real(nk-k) +   &
          (real(nj-j)*xt(5) + real(j-1)*xt(7)) * real(k-1)  )   &
       /  (real(nj-1)*real(nk-1))
    yf2=( (real(nj-j)*yt(1) + real(j-1)*yt(3)) * real(nk-k) +   &
          (real(nj-j)*yt(5) + real(j-1)*yt(7)) * real(k-1)  )   &
       /  (real(nj-1)*real(nk-1))
    zf2=( (real(nj-j)*zt(1) + real(j-1)*zt(3)) * real(nk-k) +   &
          (real(nj-j)*zt(5) + real(j-1)*zt(7)) * real(k-1)  )   &
       /  (real(nj-1)*real(nk-1))
  else
    xf2 = grid % xn(n2)
    yf2 = grid % yn(n2)
    zf2 = grid % zn(n2)
  end if

  ! Face IV
  if(grid % xn(n4) .gt. TERA) then
    xf4=( (real(nj-j)*xt(2) + real(j-1)*xt(4)) * real(nk-k) +   &
          (real(nj-j)*xt(6) + real(j-1)*xt(8)) * real(k-1)  )   &
       /  (real(nj-1)*real(nk-1))
    yf4=( (real(nj-j)*yt(2) + real(j-1)*yt(4)) * real(nk-k) +   &
          (real(nj-j)*yt(6) + real(j-1)*yt(8)) * real(k-1)  )   &
       /  (real(nj-1)*real(nk-1))
    zf4=( (real(nj-j)*zt(2) + real(j-1)*zt(4)) * real(nk-k) +   &
          (real(nj-j)*zt(6) + real(j-1)*zt(8)) * real(k-1)  )   &
       /  (real(nj-1)*real(nk-1))
  else
    xf4 = grid % xn(n4)
    yf4 = grid % yn(n4)
    zf4 = grid % zn(n4)
  end if

  if( grid % xn(n) .gt. TERA ) then
    grid % xn(n) = ( xf1*real(nk-k) + xf6*real(k-1) ) * wx16 / real(nk-1)  &
                 + ( xf2*real(ni-i) + xf4*real(i-1) ) * wx24 / real(ni-1)  &
                 + ( xf3*real(nj-j) + xf5*real(j-1) ) * wx35 / real(nj-1)

    grid % yn(n) = ( yf1*real(nk-k) + yf6*real(k-1) ) * wy16 / real(nk-1)  &
                 + ( yf2*real(ni-i) + yf4*real(i-1) ) * wy24 / real(ni-1)  &
                 + ( yf3*real(nj-j) + yf5*real(j-1) ) * wy35 / real(nj-1)

    grid % zn(n) = ( zf1*real(nk-k) + zf6*real(k-1) ) * wz16 / real(nk-1)  &
                 + ( zf2*real(ni-i) + zf4*real(i-1) ) * wz24 / real(ni-1)  &
                 + ( zf3*real(nj-j) + zf5*real(j-1) ) * wz35 / real(nj-1)
  end if

  end subroutine
