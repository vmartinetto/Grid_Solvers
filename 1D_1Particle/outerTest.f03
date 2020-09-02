include "sparseMod.f03"
      Program outerTest
      Use sparseMod
!
!     Test file to check the cross product functions of the sparse
!     linear algebra module
!
      Implicit None
      Real, Dimension(3):: vector1, vector2
      Real, Dimension(3,3):: v1xv2, v2xv1
!
!     set the vectors 1 and 2
!
      vector1 = (/1,2,3/)
      vector2 = (/2,3,4/)
!
!     do some outer products
!
      v1xv2 = outer_product(vector1,vector2)
      v2xv1 = outer_product(vector2,vector1)
!
!     print the matricies
!
      Write(*,*) v1xv2
      Write(*,*)
      Write(*,*) v2xv1
!
      End Program outerTest
