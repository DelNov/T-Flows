!==============================================================================!
  module Side_Mod
!------------------------------------------------------------------------------!
!>  This module defines the Side_Type, primarily in Front_Mod and Surf_Mod.
!>  Side_Type includes attributes and functionalities tailored to the specific
!>  needs of these modules in relationship with meshing algorithms they use.
!>  It is interesting to note that Side_Mod is a direct descendant of a module
!>  with the same name which existed in EasyMesh, and the algorithms from
!>  EasyMesh have all been transferred to Surf_Mod, in a one-to-one fashion.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Side type   !
  !---------------!
  !> Encapsulates data needed for the usage of the Side_Type objects in
  !> Front_Mod and Surf_Mod.  In Surf_Mod, in particular, Side_Type's data
  !> members play a fundamental role in surface mesh relaxation by edge
  !> swapping algorithm, leading to improved quality of surface meshes.
  type Side_Type

    integer :: ei        !! element undefined
    integer :: ea        !! element to the left
    integer :: eb        !! element to the right
    integer :: a         !! vertex from ea, which is neither c nor d
    integer :: b         !! vertex from eb, which is neither c nor d
    integer :: c, d      !! vertex defining the side (side goes from c to d)
    real    :: length    !! lengft of the side
    logical :: boundary  !! boundary tag

  end type

  end module
