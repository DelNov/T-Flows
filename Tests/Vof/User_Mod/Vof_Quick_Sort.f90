recursive subroutine Vof_Quick_Sort(a, coords, dir, first, last)
  implicit none
  real    :: coords(:,:), x, t(3)
  integer :: a(:), first, last, ti, dir
  integer :: i, j

  x = coords( (first+last) / 2, dir )
  i = first
  j = last
  do
     do while (coords(i, dir) < x)
        i=i+1
     end do
     do while (x < coords(j, dir))
        j=j-1
     end do
     if (i >= j) exit
     ti = a(i);    a(i)   = a(j);    a(j)   = ti
     t  = coords(i,:);  coords(i,:) = coords(j,:);  coords(j,:) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call Vof_Quick_Sort(a, coords, dir, first, i-1)
  if (j+1 < last)  call Vof_Quick_Sort(a, coords, dir, j+1, last)
end subroutine
