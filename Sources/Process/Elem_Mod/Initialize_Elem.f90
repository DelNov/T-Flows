!==============================================================================!
  subroutine Initialize_Elem(Elem)
!------------------------------------------------------------------------------!
!   Initializes element, as the name clearly implies                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Elem_Type) :: Elem
!==============================================================================!

  Elem % nne  = 0
  Elem % nv   = 0
  Elem % ns   = 0
  Elem % v(:) = 0
  Elem % e(:) = 0
  Elem % s(:) = 0
  Elem % cell = 0
  Elem % face = 0
  Elem % nx   = 0.0
  Elem % ny   = 0.0
  Elem % nz   = 0.0
  Elem % xc   = 0.0
  Elem % yc   = 0.0
  Elem % zc   = 0.0
  Elem % xe   = 0.0
  Elem % ye   = 0.0
  Elem % ze   = 0.0
  Elem % sx   = 0.0
  Elem % sy   = 0.0
  Elem % sz   = 0.0
  Elem % area = 0.0
  Elem % curv = 0.0

  end subroutine
