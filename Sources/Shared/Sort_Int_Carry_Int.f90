!==============================================================================!
  subroutine Sort_Int_Carry_Int(x, y, n, kflag)
!------------------------------------------------------------------------------!
  implicit  none    
!------------------------------------------------------------------------------!
!   Sorts int. array x and makes the same changes in int. array y.             !
!   It was downladed from nist and then slightly modified.                     !
!------------------------------------------------------------------------------!
  integer:: n
  integer, dimension (1:21) :: il, iu
  integer:: x(n), y(n), t, tt, ty, tty, kflag, i, j, kk, m, nn, l, ij, k
  real   :: r
!==============================================================================!

      nn = n
      if (nn.ge.1) goto 10
      print *, 'Sort_Int_Carry_Int: the number of values to be ' // &
               'sorted was not positive.'
      return
   10 kk = iabs(kflag)
      if ((kk.eq.1).or.(kk.eq.2)) goto 15
      print *, 'Sort_Int_Carry_Int: the sort control parameter, ' // &
               'k, was not 2, 1, -1, or -2.'
      return

!-----------------------------------------------------!
!   Alter array x to get decreasing order if needed   !
!-----------------------------------------------------!
   15 if (kflag.ge.1) goto 30
      do 20 i=1,nn
      x(i) = -x(i)
   20 continue
   30 if(kk .eq. 1) goto 100
      if(kk .eq. 2) goto 200 

!-----------------!
!   Sort x only   !
!-----------------!
  100 continue
      m=1
      i=1
      j=nn
      r=.375
  110 if (i .eq. j) goto 155
      if (r .gt. .5898437) goto 120
      r=r+3.90625e-2
      goto 125
  120 r=r-.21875
  125 k=i
!                                  select a central element of the
!                                  array and save it in location t
      ij = i + ifix (float (j-i) * r)
      t=x(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) goto 130
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
  130 l=j
!                                  if last element of array is less than
!                                  t, interchange with t
      if (x(j) .ge. t) goto 140
      x(ij)=x(j)
      x(j)=t
      t=x(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) goto 140
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
      goto 140
  135 tt=x(l)
      x(l)=x(k)
      x(k)=tt
!                                  find an element in the second half of
!                                  the array which is smaller than t
  140 l=l-1
      if (x(l) .gt. t) goto 140
!                                  find an element in the first half of
!                                  the array which is greater than t
  145 k=k+1
      if (x(k) .lt. t) goto 145
!                                  interchange these elements
      if (k .le. l) goto 135
!                                  save upper and lower subscripts of
!                                  the array yet to be sorted
      if (l-i .le. j-k) goto 150
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      goto 160
  150 il(m)=k
      iu(m)=j
      j=l
      m=m+1
      goto 160
!                                  begin again on another portion of
!                                  the unsorted array
  155 m=m-1
      if (m .eq. 0) goto 300
      i=il(m)
      j=iu(m)
  160 if (j-i .ge. 1) goto 125
      if (i .eq. 1) goto 110
      i=i-1
  165 i=i+1
      if (i .eq. j) goto 155
      t=x(i+1)
      if (x(i) .le. t) goto 165
      k=i
  170 x(k+1)=x(k)
      k=k-1
      if (t .lt. x(k)) goto 170
      x(k+1)=t
      goto 165

!------------------------------!
!   Sort x and carry y along   !
!------------------------------!
  200 continue
      m=1
      i=1
      j=nn
      r=.375
  210 if (i .eq. j) goto 255
      if (r .gt. .5898437) goto 220
      r=r+3.90625e-2
      goto 225
  220 r=r-.21875
  225 k=i
!                                  select a central element of the
!                                  array and save it in location t
      ij = i + ifix (float (j-i) *r)
      t=x(ij)
      ty= y(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) goto 230
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
       y(ij)= y(i)
       y(i)=ty
      ty= y(ij)
  230 l=j
!                                  if last element of array is less than
!                                  t, interchange with t
      if (x(j) .ge. t) goto 240
      x(ij)=x(j)
      x(j)=t
      t=x(ij)
       y(ij)= y(j)
       y(j)=ty
      ty= y(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (x(i) .le. t) goto 240
      x(ij)=x(i)
      x(i)=t
      t=x(ij)
       y(ij)= y(i)
       y(i)=ty
      ty= y(ij)
      goto 240
  235 tt=x(l)
      x(l)=x(k)
      x(k)=tt
      tty= y(l)
       y(l)= y(k)
       y(k)=tty
!                                  find an element in the second half of
!                                  the array which is smaller than t
  240 l=l-1
      if (x(l) .gt. t) goto 240
!                                  find an element in the first half of
!                                  the array which is greater than t
  245 k=k+1
      if (x(k) .lt. t) goto 245
!                                  interchange these elements
      if (k .le. l) goto 235
!                                  save upper and lower subscripts of
!                                  the array yet to be sorted
      if (l-i .le. j-k) goto 250
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      goto 260
  250 il(m)=k
      iu(m)=j
      j=l
      m=m+1
      goto 260
!                                  begin again on another portion of
!                                  the unsorted array
  255 m=m-1
      if (m .eq. 0) goto 300
      i=il(m)
      j=iu(m)
  260 if (j-i .ge. 1) goto 225
      if (i .eq. 1) goto 210
      i=i-1
  265 i=i+1
      if (i .eq. j) goto 255
      t=x(i+1)
      ty= y(i+1)
      if (x(i) .le. t) goto 265
      k=i
  270 x(k+1)=x(k)
       y(k+1)= y(k)
      k=k-1
      if (t .lt. x(k)) goto 270
      x(k+1)=t
       y(k+1)=ty
      goto 265

!--------------!
!   Clean up   !
!--------------!
  300 if (kflag.ge.1) return
      do 310 i=1,nn
      x(i) = -x(i)
  310 continue

      end subroutine
